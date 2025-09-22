#include "text_io.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <vector>

#include <zlib.h>
#include <zstd.h>

namespace {

class PlainLineReader : public LineReader {
public:
    explicit PlainLineReader(const std::string& path) {
        stream_.open(path);
        if (!stream_.is_open()) {
            throw std::runtime_error("Failed to open file " + path);
        }
    }

    bool getline(std::string& line) override {
        return static_cast<bool>(std::getline(stream_, line));
    }

private:
    std::ifstream stream_;
};

class GzipLineReader : public LineReader {
public:
    explicit GzipLineReader(const std::string& path) {
        file_ = gzopen(path.c_str(), "rb");
        if (file_ == Z_NULL) {
            throw std::runtime_error("Failed to open gzip file " + path);
        }
        buffer_.resize(kBufferSize);
    }

    ~GzipLineReader() override {
        if (file_ != Z_NULL) {
            gzclose(file_);
        }
    }

    bool getline(std::string& line) override {
        line.clear();
        if (file_ == Z_NULL) {
            return false;
        }

        while (true) {
            char* result = gzgets(file_, buffer_.data(), static_cast<int>(buffer_.size()));
            if (result == Z_NULL) {
                if (gzeof(file_) && !line.empty()) {
                    return true;
                }
                return false;
            }
            size_t len = std::strlen(result);
            bool has_newline = len > 0 && result[len - 1] == '\n';
            if (has_newline) {
                len -= 1;
            }
            line.append(result, len);
            if (has_newline) {
                return true;
            }
        }
    }

private:
    static constexpr size_t kBufferSize = 1 << 15;
    gzFile file_ = Z_NULL;
    std::vector<char> buffer_;
};

class ZstdLineReader : public LineReader {
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

    ~ZstdLineReader() override {
        if (dctx_) {
            ZSTD_freeDCtx(dctx_);
        }
        if (file_) {
            std::fclose(file_);
        }
    }

    bool getline(std::string& line) override {
        while (true) {
            const auto newline_pos = pending_.find('\n');
            if (newline_pos != std::string::npos) {
                line = pending_.substr(0, newline_pos);
                pending_.erase(0, newline_pos + 1);
                return true;
            }
            if (stream_finished_) {
                if (pending_.empty()) {
                    return false;
                }
                line.swap(pending_);
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

class PlainWriterImpl {
public:
    explicit PlainWriterImpl(const std::string& path) {
        stream_.open(path, std::ios::out | std::ios::binary);
        if (!stream_.is_open()) {
            throw std::runtime_error("Failed to open file for writing: " + path);
        }
    }

    void write(const char* data, size_t size) {
        stream_.write(data, static_cast<std::streamsize>(size));
        if (!stream_) {
            throw std::runtime_error("Failed to write to file");
        }
    }

    void flush() {
        stream_.flush();
    }

private:
    std::ofstream stream_;
};

class GzipWriterImpl {
public:
    explicit GzipWriterImpl(const std::string& path) {
        file_ = gzopen(path.c_str(), "wb");
        if (file_ == Z_NULL) {
            throw std::runtime_error("Failed to open gzip file for writing: " + path);
        }
    }

    ~GzipWriterImpl() {
        if (file_ != Z_NULL) {
            gzclose(file_);
        }
    }

    void write(const char* data, size_t size) {
        if (gzwrite(file_, data, static_cast<unsigned int>(size)) != static_cast<int>(size)) {
            throw std::runtime_error("Failed to write to gzip file");
        }
    }

    void flush() {
        gzflush(file_, Z_SYNC_FLUSH);
    }

private:
    gzFile file_ = Z_NULL;
};

class ZstdWriterImpl {
public:
    explicit ZstdWriterImpl(const std::string& path) {
        file_ = std::fopen(path.c_str(), "wb");
        if (!file_) {
            throw std::runtime_error("Failed to open zstd file for writing: " + path);
        }
        cctx_ = ZSTD_createCCtx();
        if (!cctx_) {
            std::fclose(file_);
            throw std::runtime_error("Failed to create zstd compression context");
        }
        out_buffer_.resize(kBufferSize);
    }

    ~ZstdWriterImpl() {
        finish();
        if (cctx_) {
            ZSTD_freeCCtx(cctx_);
        }
        if (file_) {
            std::fclose(file_);
        }
    }

    void write(const char* data, size_t size) {
        if (finished_) {
            throw std::runtime_error("Cannot write after finalizing ZstdWriter");
        }
        ZSTD_inBuffer input{data, size, 0};
        while (input.pos < input.size) {
            ZSTD_outBuffer output{out_buffer_.data(), out_buffer_.size(), 0};
            size_t ret = ZSTD_compressStream2(cctx_, &output, &input, ZSTD_e_continue);
            if (ZSTD_isError(ret)) {
                throw std::runtime_error(std::string("ZSTD compression error: ") + ZSTD_getErrorName(ret));
            }
            flush_output(output);
        }
    }

    void finish() {
        if (finished_) {
            return;
        }
        ZSTD_inBuffer input{nullptr, 0, 0};
        while (true) {
            ZSTD_outBuffer output{out_buffer_.data(), out_buffer_.size(), 0};
            size_t ret = ZSTD_compressStream2(cctx_, &output, &input, ZSTD_e_end);
            if (ZSTD_isError(ret)) {
                throw std::runtime_error(std::string("ZSTD finalize error: ") + ZSTD_getErrorName(ret));
            }
            flush_output(output);
            if (ret == 0) {
                break;
            }
        }
        finished_ = true;
    }

    void flush() {
        if (finished_) {
            return;
        }
        std::fflush(file_);
    }

private:
    void flush_output(ZSTD_outBuffer& output) {
        if (output.pos > 0) {
            size_t written = std::fwrite(out_buffer_.data(), 1, output.pos, file_);
            if (written != output.pos) {
                throw std::runtime_error("Failed to write compressed data to file");
            }
        }
    }

    static constexpr size_t kBufferSize = 1 << 16;
    std::FILE* file_ = nullptr;
    ZSTD_CCtx* cctx_ = nullptr;
    std::vector<char> out_buffer_;
    bool finished_ = false;
};

} // namespace

struct TextWriter::Impl {
    Impl() = default;

    std::unique_ptr<PlainWriterImpl> plain;
    std::unique_ptr<GzipWriterImpl> gzip;
    std::unique_ptr<ZstdWriterImpl> zstd;
    CompressionType type = CompressionType::PLAIN;

    void write(const char* data, size_t size) {
        switch (type) {
            case CompressionType::PLAIN:
                plain->write(data, size);
                break;
            case CompressionType::GZIP:
                gzip->write(data, size);
                break;
            case CompressionType::ZSTD:
                zstd->write(data, size);
                break;
            case CompressionType::AUTO:
                throw std::logic_error("CompressionType::AUTO is invalid for TextWriter");
        }
    }

    void flush() {
        switch (type) {
            case CompressionType::PLAIN:
                plain->flush();
                break;
            case CompressionType::GZIP:
                gzip->flush();
                break;
            case CompressionType::ZSTD:
                zstd->flush();
                break;
            case CompressionType::AUTO:
                break;
        }
    }
};

TextWriter::TextWriter() : impl_(new Impl()) {}

TextWriter::~TextWriter() = default;

std::unique_ptr<TextWriter> TextWriter::create(const std::string& path, CompressionType type) {
    auto writer = std::unique_ptr<TextWriter>(new TextWriter());
    writer->impl_->type = type;
    switch (type) {
        case CompressionType::PLAIN:
            writer->impl_->plain.reset(new PlainWriterImpl(path));
            break;
        case CompressionType::GZIP:
            writer->impl_->gzip.reset(new GzipWriterImpl(path));
            break;
        case CompressionType::ZSTD:
            writer->impl_->zstd.reset(new ZstdWriterImpl(path));
            break;
        case CompressionType::AUTO:
            throw std::logic_error("CompressionType::AUTO is invalid for TextWriter");
    }
    return writer;
}

void TextWriter::write(const std::string& data) {
    impl_->write(data.data(), data.size());
}

void TextWriter::write(const char* data, size_t size) {
    impl_->write(data, size);
}

void TextWriter::flush() {
    impl_->flush();
}

std::unique_ptr<LineReader> make_line_reader(const std::string& path, CompressionType type) {
    switch (type) {
        case CompressionType::PLAIN:
            return std::unique_ptr<LineReader>(new PlainLineReader(path));
        case CompressionType::GZIP:
            return std::unique_ptr<LineReader>(new GzipLineReader(path));
        case CompressionType::ZSTD:
            return std::unique_ptr<LineReader>(new ZstdLineReader(path));
        case CompressionType::AUTO:
            throw std::logic_error("CompressionType::AUTO is invalid for make_line_reader");
    }
    throw std::logic_error("Unhandled compression type");
}
