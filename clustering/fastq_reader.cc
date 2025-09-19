#include "fastq_reader.h"

#include "kseq.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <vector>

#include <zlib.h>

#include <zstd.h>

namespace {

KSEQ_INIT(gzFile, gzread)

[[maybe_unused]] auto* kseq_rewind_ptr = &kseq_rewind;

std::string toLower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

bool hasSuffix(const std::string& value, const std::string& suffix) {
    if (suffix.size() > value.size()) {
        return false;
    }
    return std::equal(suffix.rbegin(), suffix.rend(), value.rbegin(), [](char a, char b) {
        return std::tolower(static_cast<unsigned char>(a)) == std::tolower(static_cast<unsigned char>(b));
    });
}

class ZstdLineReader {
public:
    explicit ZstdLineReader(const std::string& path) {
        file_ = std::fopen(path.c_str(), "rb");
        if (!file_) {
            throw std::runtime_error("Failed to open zstd file: " + path);
        }
        dctx_ = ZSTD_createDCtx();
        if (!dctx_) {
            std::fclose(file_);
            throw std::runtime_error("Failed to create zstd decompression context");
        }
        in_buffer_.resize(kBufferSize);
        out_buffer_.resize(kBufferSize);
        input_.src = in_buffer_.data();
        input_.size = 0;
        input_.pos = 0;
    }

    ~ZstdLineReader() {
        if (dctx_) {
            ZSTD_freeDCtx(dctx_);
        }
        if (file_) {
            std::fclose(file_);
        }
    }

    bool getline(std::string& line) {
        while (true) {
            const auto newline_pos = pending_.find('\n');
            if (newline_pos != std::string::npos) {
                line = pending_.substr(0, newline_pos);
                pending_.erase(0, newline_pos + 1);
                if (!line.empty() && line.back() == '\r') {
                    line.pop_back();
                }
                return true;
            }
            if (stream_finished_) {
                if (pending_.empty()) {
                    return false;
                }
                line.swap(pending_);
                if (!line.empty() && line.back() == '\r') {
                    line.pop_back();
                }
                pending_.clear();
                return true;
            }
            pump();
        }
    }

private:
    void pump() {
        if (stream_finished_) {
            return;
        }
        if (input_.pos == input_.size && !input_depleted_) {
            const size_t bytes_read = std::fread(in_buffer_.data(), 1, in_buffer_.size(), file_);
            input_.src = in_buffer_.data();
            input_.size = bytes_read;
            input_.pos = 0;
            if (bytes_read == 0) {
                input_depleted_ = true;
            }
        }

        ZSTD_outBuffer output;
        output.dst = out_buffer_.data();
        output.size = out_buffer_.size();
        output.pos = 0;

        const size_t ret = ZSTD_decompressStream(dctx_, &output, &input_);
        if (ZSTD_isError(ret)) {
            throw std::runtime_error(std::string("ZSTD decompression error: ") + ZSTD_getErrorName(ret));
        }

        if (output.pos > 0) {
            pending_.append(out_buffer_.data(), output.pos);
        }

        if (ret == 0 && input_depleted_ && input_.pos == input_.size) {
            stream_finished_ = true;
        }
        if (output.pos == 0 && input_depleted_ && input_.pos == input_.size) {
            stream_finished_ = true;
        }
    }

    static constexpr size_t kBufferSize = 1 << 16;

    std::FILE* file_ = nullptr;
    ZSTD_DCtx* dctx_ = nullptr;
    std::vector<char> in_buffer_;
    std::vector<char> out_buffer_;
    ZSTD_inBuffer input_{};
    bool input_depleted_ = false;
    bool stream_finished_ = false;
    std::string pending_;
};

} // namespace

CompressionType detectCompression(const std::string& path) {
    const std::string lower = toLower(path);
    if (hasSuffix(lower, ".zst") || hasSuffix(lower, ".zstd")) {
        return CompressionType::ZSTD;
    }
    if (hasSuffix(lower, ".gz") || hasSuffix(lower, ".gzip")) {
        return CompressionType::GZIP;
    }
    return CompressionType::PLAIN;
}

std::string compressionTypeName(CompressionType type) {
    switch (type) {
        case CompressionType::AUTO:
            return "auto";
        case CompressionType::PLAIN:
            return "plain";
        case CompressionType::GZIP:
            return "gzip";
        case CompressionType::ZSTD:
            return "zstd";
        default:
            return "unknown";
    }
}

CompressionType resolveCompression(CompressionType requested, const std::string& path) {
    if (requested == CompressionType::AUTO) {
        return detectCompression(path);
    }
    return requested;
}

class FastqReader::Impl {
public:
    Impl(const std::string& path, CompressionType compression)
        : compression_(resolveCompression(compression, path)), path_(path) {
        switch (compression_) {
            case CompressionType::PLAIN:
                stream_.open(path);
                if (!stream_.is_open()) {
                    throw std::runtime_error("Failed to open FASTQ file: " + path);
                }
                break;
            case CompressionType::GZIP:
                gz_stream_ = gzopen(path.c_str(), "rb");
                if (gz_stream_ == Z_NULL) {
                    throw std::runtime_error("Failed to open gzipped FASTQ file: " + path);
                }
                kseq_reader_ = kseq_init(gz_stream_);
                if (kseq_reader_ == nullptr) {
                    gzclose(gz_stream_);
                    throw std::runtime_error("Failed to create kseq reader for file: " + path);
                }
                break;
            case CompressionType::ZSTD:
                zstd_reader_.reset(new ZstdLineReader(path));
                break;
            case CompressionType::AUTO:
                throw std::logic_error("CompressionType::AUTO should not reach FastqReader implementation");
        }
    }

    ~Impl() {
        if (kseq_reader_ != nullptr) {
            kseq_destroy(kseq_reader_);
        }
        if (gz_stream_ != Z_NULL) {
            gzclose(gz_stream_);
        }
        if (stream_.is_open()) {
            stream_.close();
        }
    }

    bool readRecord(FastqRecord& record) {
        switch (compression_) {
            case CompressionType::PLAIN:
                return readPlain(record);
            case CompressionType::GZIP:
                return readGzip(record);
            case CompressionType::ZSTD:
                return readZstd(record);
            case CompressionType::AUTO:
                return false;
        }
        return false;
    }

private:
    bool readPlain(FastqRecord& record) {
        if (!stream_.good()) {
            return false;
        }
        if (!std::getline(stream_, record.name)) {
            return false;
        }
        if (!std::getline(stream_, record.sequence)) {
            return false;
        }
        if (!std::getline(stream_, record.plus)) {
            return false;
        }
        if (!std::getline(stream_, record.quality)) {
            return false;
        }
        return true;
    }

    bool readGzip(FastqRecord& record) {
        if (kseq_reader_ == nullptr) {
            return false;
        }
        const int status = kseq_read(kseq_reader_);
        if (status < 0) {
            return false;
        }
        record.name.assign("@");
        record.name.append(kseq_reader_->name.s, kseq_reader_->name.l);
        record.sequence.assign(kseq_reader_->seq.s, kseq_reader_->seq.l);
        record.plus = "+";
        if (kseq_reader_->qual.l > 0) {
            record.quality.assign(kseq_reader_->qual.s, kseq_reader_->qual.l);
        } else {
            record.quality.clear();
        }
        return true;
    }

    bool readZstd(FastqRecord& record) {
        if (!zstd_reader_) {
            return false;
        }
        if (!zstd_reader_->getline(record.name)) {
            return false;
        }
        if (!zstd_reader_->getline(record.sequence)) {
            return false;
        }
        if (!zstd_reader_->getline(record.plus)) {
            return false;
        }
        if (!zstd_reader_->getline(record.quality)) {
            return false;
        }
        return true;
    }

    CompressionType compression_;
    std::string path_;
    std::ifstream stream_;
    gzFile gz_stream_ = Z_NULL;
    kseq_t* kseq_reader_ = nullptr;
    std::unique_ptr<ZstdLineReader> zstd_reader_;
};

FastqReader::FastqReader(const std::string& path, CompressionType compression)
    : impl_(new Impl(path, compression)) {}

FastqReader::~FastqReader() {
    delete impl_;
}

bool FastqReader::readRecord(FastqRecord& record) {
    return impl_->readRecord(record);
}
