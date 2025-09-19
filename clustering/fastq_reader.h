#ifndef CALIB_FASTQ_READER_H
#define CALIB_FASTQ_READER_H

#include <string>

enum class CompressionType {
    AUTO = 0,
    PLAIN,
    GZIP,
    ZSTD
};

CompressionType detectCompression(const std::string& path);
std::string compressionTypeName(CompressionType type);

struct FastqRecord {
    std::string name;
    std::string sequence;
    std::string plus;
    std::string quality;
};

class FastqReader {
public:
    FastqReader(const std::string& path, CompressionType compression);
    ~FastqReader();

    FastqReader(const FastqReader&) = delete;
    FastqReader& operator=(const FastqReader&) = delete;

    bool readRecord(FastqRecord& record);

private:
    class Impl;
    Impl* impl_;
};

CompressionType resolveCompression(CompressionType requested, const std::string& path);

#endif // CALIB_FASTQ_READER_H
