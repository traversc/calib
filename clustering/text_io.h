#ifndef CALIB_TEXT_IO_H
#define CALIB_TEXT_IO_H

#include <memory>
#include <string>

#include "fastq_reader.h"

class LineReader {
public:
    virtual ~LineReader() = default;
    virtual bool getline(std::string& line) = 0;
};

std::unique_ptr<LineReader> make_line_reader(const std::string& path, CompressionType type);

class TextWriter {
public:
    static std::unique_ptr<TextWriter> create(const std::string& path, CompressionType type);

    TextWriter(const TextWriter&) = delete;
    TextWriter& operator=(const TextWriter&) = delete;

    void write(const std::string& data);
    void write(const char* data, size_t size);
    void flush();

    ~TextWriter();

private:
    TextWriter();
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

#endif // CALIB_TEXT_IO_H
