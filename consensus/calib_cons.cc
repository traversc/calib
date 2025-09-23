#include "spoa/spoa.hpp"
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <iterator>
#include <fstream>
#include <iostream>
#include <thread>
#include <functional>
#include <locale>
#include <exception>
#include <algorithm>
#include <cctype>
#include <array>
#include <cmath>
#include <vector>

#include "../clustering/fastq_reader.h"
#include "../clustering/text_io.h"

#define ALIGNMENT_TYPE_GLOBAL 1
#define SCORE_M 5
#define SCORE_X -3
#define SCORE_G -9
#define ASCII_SIZE 128
#define MSA_MAJORITY 0.5

int thread_count = 4;
int min_reads_per_cluster = 2;
int max_reads_per_cluster = 1000;
std::string cluster_filename = "";
std::vector<std::string> fastq_filenames;
std::vector<std::string> output_filenames;
CompressionType fastq_compression = CompressionType::AUTO;
CompressionType cluster_input_compression = CompressionType::AUTO;
CompressionType output_compression = CompressionType::PLAIN;
typedef uint32_t cluster_id_t;
typedef uint32_t read_id_t;


void print_help(){
    std::cout << "Calib Consensus: Generating consensus sequence from Calib clusters." << "\n";
    std::cout << "Usage: calib_cons [--PARAMETER VALUE]" << "\n";
    std::cout << "Example 1: calib_cons -t 8 -c input.cluster -q 1.fastq 2.fastq -o 1.out 2.out" << "\n";
    std::cout << "Example 2: calib_cons -q 1.fastq -q 2.fastq -o 1.out 2.out -c input.cluster" << "\n";
    std::cout << "Calib's paramters arguments:" << "\n";
    std::cout << "  -q  --fastq                    (type: space separated string list;\n";
    std::cout << "                                    REQUIRED paramter;\n";
    std::cout << "                                    can be set multiple times like in Example 2)\n";
    std::cout << "  -o  --output-prefix            (type: space separated string list;\n";
    std::cout << "                                    REQUIRED paramter;\n";
    std::cout << "                                    can be set multiple times like in Example 2;\n";
    std::cout << "                                    must be same size as fastq list)\n";
    std::cout << "  -c  --cluster                  (string;\n";
    std::cout << "                                    REQUIRED paramter)\n";
    std::cout << "  -t  --threads                  (positive integer;\n";
    std::cout << "                                    default: 4)\n";
    std::cout << "  -m  --min-reads-per-cluster    (positive integer;\n";
    std::cout << "                                    default: 2)\n";
    std::cout << "  -x  --max-reads-per-cluster    (positive integer;\n";
    std::cout << "                                    default: 1000)\n";
    std::cout << "      --input-format <auto|plain|gzip|zstd>      (type: string; default: auto)\n";
    std::cout << "      --cluster-input-format <auto|plain|gzip|zstd> (type: string; default: auto)\n";
    std::cout << "      --output-format <plain|gzip|zstd>          (type: string; default: plain)\n";
    std::cout << "  -h  --help\n";
}

void parse_flags(int argc, char *argv[]){
    bool input_format_specified = false;
    bool output_format_specified = false;
    auto parse_format_value = [&](const std::string& value, const std::string& flag, bool allow_auto) {
        std::string lower = value;
        std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
        if (allow_auto && lower == "auto") {
            return CompressionType::AUTO;
        }
        if (lower == "plain") {
            return CompressionType::PLAIN;
        }
        if (lower == "gzip") {
            return CompressionType::GZIP;
        }
        if (lower == "zstd" || lower == "zstandard") {
            return CompressionType::ZSTD;
        }
        std::cout << "Unknown format value '" << value << "' for " << flag << "\n";
        print_help();
        exit(-1);
        return CompressionType::AUTO;
    };
    for (int i = 1; i < argc; i++) {
        std::string current_param(argv[i]);
        if (current_param == "-h" || current_param == "--help") {
            print_help();
            exit(0);
        }
        if (current_param == "-c" || current_param == "--cluster") {
            i++;
            cluster_filename = std::string(argv[i]);
            continue;
        }
        if (current_param == "-t" || current_param == "--threads") {
            i++;
            thread_count = atoi(argv[i]);
            continue;
        }
        if (current_param == "-m" || current_param == "--min-reads-per-cluster") {
            i++;
            min_reads_per_cluster = atoi(argv[i]);
            continue;
        }
        if (current_param == "-x" || current_param == "--max-reads-per-cluster") {
            i++;
            max_reads_per_cluster = atoi(argv[i]);
            continue;
        }
        if (current_param == "-q" || current_param == "--fastq") {
            i++;
            for (; i < argc; i++) {
                std::string fastq_filename(argv[i]);
                if (fastq_filename[0] == '-') {
                    i--;
                    break;
                }
                fastq_filenames.push_back(fastq_filename);
            }
            continue;
        }
        if (current_param == "-o" || current_param == "--output-prefix") {
            i++;
            for (; i < argc; i++) {
                std::string output_filename(argv[i]);
                if (output_filename[0] == '-') {
                    i--;
                    break;
                }
                output_filenames.push_back(output_filename);
            }
            continue;
        }
        if (current_param == "--input-format") {
            if (i + 1 >= argc) {
                std::cout << "Missing value for --input-format" << '\n';
                print_help();
                exit(-1);
            }
            CompressionType requested = parse_format_value(std::string(argv[i+1]), current_param, true);
            if (input_format_specified && requested != fastq_compression) {
                std::cout << "Conflicting input format. Already set to "
                          << compressionTypeName(fastq_compression)
                          << ", cannot change to " << compressionTypeName(requested)
                          << " via " << current_param << "\n";
                print_help();
                exit(-1);
            }
            fastq_compression = requested;
            input_format_specified = true;
            i++;
            continue;
        }
        if (current_param == "--cluster-input-format") {
            if (i + 1 >= argc) {
                std::cout << "Missing value for --cluster-input-format" << '\n';
                print_help();
                exit(-1);
            }
            cluster_input_compression = parse_format_value(std::string(argv[i+1]), current_param, true);
            i++;
            continue;
        }
        if (current_param == "--output-format") {
            if (i + 1 >= argc) {
                std::cout << "Missing value for --output-format" << '\n';
                print_help();
                exit(-1);
            }
            CompressionType requested = parse_format_value(std::string(argv[i+1]), current_param, false);
            if (output_format_specified && requested != output_compression) {
                std::cout << "Conflicting output format. Already set to "
                          << compressionTypeName(output_compression)
                          << ", cannot change to " << compressionTypeName(requested)
                          << " via " << current_param << "\n";
                print_help();
                exit(-1);
            }
            output_compression = requested;
            output_format_specified = true;
            i++;
            continue;
        }
        std::cout << "Unrecognized parameter, " << argv[i] << ", was passed.\n";
        print_help();
        exit(-1);
    }
    if (min_reads_per_cluster < 0 ) {
        std::cout << "Minimum reads per cluster ("<<min_reads_per_cluster<<") must be positive value." << '\n';
        print_help();
        exit(-1);
    }
    if (max_reads_per_cluster < 0 ) {
        std::cout << "Minimum reads per cluster ("<<max_reads_per_cluster<<") must be positive value." << '\n';
        print_help();
        exit(-1);
    }
    if (thread_count < 0 || thread_count > 16) {
        std::cout << "Thread count must be between 1 and 16." << '\n';
        print_help();
        exit(-1);
    }
    if (fastq_filenames.size() != output_filenames.size()) {
        std::cout << "Number of fastq files ("<< fastq_filenames.size() <<") must be equal number of output files ("<< output_filenames.size() <<")\n";
        print_help();
        exit(-1);
    }
    if (cluster_filename == "") {
        std::cout << "No cluster filename was passed.\n";
        print_help();
        exit(-1);
    }

}

char consensus_quality_from_agreement(double agreement_fraction) {
    if (!std::isfinite(agreement_fraction) || agreement_fraction <= 0.0) {
        return static_cast<char>('!' + 0);
    }
    double clamped = std::max(0.5, std::min(1.0, agreement_fraction));
    double normalized = (clamped - 0.5) / 0.5;
    int phred_value = static_cast<int>(std::round(normalized * 93.0));
    phred_value = std::max(0, std::min(93, phred_value));
    return static_cast<char>('!' + phred_value);
}

void process_clusters(const std::vector<std::string>& read_to_sequence,
                      const std::vector<std::string>& read_to_quality,
                      const std::vector<std::vector<read_id_t> >& cluster_to_reads,
                      std::string o_filename_prefix,
                      size_t thread_id) {
    std::ofstream ofastq(o_filename_prefix + "fastq" + std::to_string(thread_id));
    std::ofstream omsa(o_filename_prefix + "msa" + std::to_string(thread_id));
    for (cluster_id_t cid = thread_id; cid < cluster_to_reads.size(); cid+=thread_count) {
        if (cluster_to_reads[cid].size() < min_reads_per_cluster) {
            continue;
        }
        if (cluster_to_reads[cid].size() > max_reads_per_cluster) {
            continue;
        }
        std::stringstream header;
        header << cid << "\t";

        const auto& read_ids = cluster_to_reads[cid];
        size_t cluster_size = read_ids.size();

        size_t reads_to_use = cluster_size;
        if (max_reads_per_cluster > 1000 && reads_to_use > 1000) {
            reads_to_use = 1000;
        }

        for (size_t idx = 0; idx < reads_to_use; ++idx) {
            header << read_ids[idx] << ";";
        }

        // If only one read in cluster, skip MSA and use the read as consensus
        if (min_reads_per_cluster == 1 && cluster_size == 1) {
            read_id_t rid = read_ids[0];
            const std::string& single_sequence = read_to_sequence[rid];
            const std::string& single_quality = read_to_quality[rid];

            ofastq << "@" << header.str() << '\n';
            ofastq << single_sequence << '\n';
            ofastq << '+' << '\n';
            ofastq << single_quality << '\n';

            omsa << "@" << header.str() << '\n';
            omsa << single_sequence << '\n';
            omsa << '+' << '\n';
            omsa << single_sequence << '\n';
        } else {
            auto alignment_engine = spoa::createAlignmentEngine(
                static_cast<spoa::AlignmentType>(ALIGNMENT_TYPE_GLOBAL), SCORE_M, SCORE_X, SCORE_G
            );
            auto graph = spoa::createGraph();
            for (size_t idx = 0; idx < reads_to_use; ++idx) {
                read_id_t rid = read_ids[idx];
                auto alignment = alignment_engine->align_sequence_with_graph(read_to_sequence[rid], graph);
                graph->add_alignment(alignment, read_to_sequence[rid]);
            }
            std::string consensus;
            std::string qual;
            std::vector<std::string> msa;
            graph->generate_multiple_sequence_alignment(msa);

            size_t profile_width = msa[0].size();
            size_t profile_height = msa.size();
            struct ColumnProfile {
                std::array<double, ASCII_SIZE> counts{};
            };
            std::vector<ColumnProfile> profile(profile_width);
            for (const auto& it : msa) {
                for (size_t col = 0; col < profile_width; col++) {
                    unsigned char symbol = static_cast<unsigned char>(it[col]);
                    unsigned char upper_symbol = symbol;
                    int upper_value = std::toupper(static_cast<unsigned char>(symbol));
                    if (upper_value >= 0 && upper_value < ASCII_SIZE) {
                        upper_symbol = static_cast<unsigned char>(upper_value);
                    }
                    profile[col].counts[upper_symbol] += 1.0;
                }
            }
            
            consensus.reserve(profile_width);
            const std::array<char, 4> bases{'A', 'C', 'G', 'T'};
            for (size_t col = 0; col < profile_width; col++) {
                if (profile_height == 0) {
                    continue;
                }
                double gap_fraction = profile[col].counts[static_cast<unsigned char>('-')] / static_cast<double>(profile_height);
                if (gap_fraction > MSA_MAJORITY) { continue; } // skip columns with majority gaps
                double best_fraction = 0.0;
                char best_base = 'N';
                for (char base : bases) {
                    double fraction = profile[col].counts[static_cast<unsigned char>(base)] / static_cast<double>(profile_height);
                    if (fraction > best_fraction) {
                        best_fraction = fraction;
                        best_base = base;
                    }
                }
                if (best_fraction > MSA_MAJORITY) {
                    consensus += best_base;
                } else {
                    consensus += 'N';
                }
                qual += consensus_quality_from_agreement(best_fraction);
            }
            ofastq << "@" << header.str() << '\n';
            ofastq << consensus << '\n';
            ofastq << '+' << '\n';
            ofastq << qual << '\n';

            omsa << "@" << header.str() << '\n';
            omsa << consensus << '\n';
            omsa << '+' << '\n';
            for (const auto& it: msa) {
                omsa << it << '\n';
            }
        }
    }
}

void run_consensus(){
    std::vector<std::vector<read_id_t> > cluster_to_reads;
    size_t read_count = 0;
    std::string line_buffer;
    CompressionType cluster_type = resolveCompression(cluster_input_compression, cluster_filename);
    auto cluster_reader = make_line_reader(cluster_filename, cluster_type);
    std::cout << "Reading cluster file: " << cluster_filename << '\n';
    while (cluster_reader->getline(line_buffer)) {
        std::stringstream ss(line_buffer);
        std::istream_iterator<std::string> begin(ss);
        std::istream_iterator<std::string> end;
        std::vector<std::string> vstrings(begin, end);
        cluster_id_t cid = atoi(vstrings[0].c_str());
        read_id_t rid = atoi(vstrings[2].c_str());

        if (cluster_to_reads.size() <= cid) {
            cluster_to_reads.resize(cid+1);
        }
        cluster_to_reads[cid].push_back(rid);
        read_count++;
    }

    //Get read sequences from each FASTQ file, and pass it for MSA and output
    for (int i = 0; i < fastq_filenames.size(); i++) {
        std::cout << "Reading fastq file: " << fastq_filenames[i] << '\n';
        std::string ifastq_filename = fastq_filenames[i];
        std::string o_filename_prefix = output_filenames[i];
        const std::string output_extension = compressionTypeExtension(output_compression);
        const std::string fastq_output_path = o_filename_prefix + "fastq" + output_extension;
        const std::string msa_output_path = o_filename_prefix + "msa" + output_extension;
        std::cout << "Writing output files: " << fastq_output_path << " and " << msa_output_path << '\n';

        std::vector<std::string> read_to_sequence(read_count);
        std::vector<std::string> read_to_quality(read_count);
        read_id_t rid = 0;
        try {
            CompressionType resolved_type = resolveCompression(fastq_compression, ifastq_filename);
            std::cout << "Detected " << compressionTypeName(resolved_type) << " compression" << '\n';
            FastqReader fastq_reader(ifastq_filename, resolved_type);
            FastqRecord record;
            while (fastq_reader.readRecord(record)) {
                if (read_to_sequence.size() <= rid) {
                    read_to_sequence.resize(rid + 1);
                    read_to_quality.resize(rid + 1);
                }
                read_to_sequence[rid] = record.sequence;
                read_to_quality[rid] = record.quality;
                rid++;
            }
        } catch (const std::exception& ex) {
            std::cout << ex.what() << '\n';
            exit(-1);
        }
        std::vector<std::thread> threads(thread_count);
        for (size_t thread_id = 0; thread_id < thread_count; thread_id++) {
            threads[thread_id] = std::thread(process_clusters,
                                            std::ref(read_to_sequence),
                                            std::ref(read_to_quality),
                                            std::ref(cluster_to_reads),
                                            o_filename_prefix,
                                            thread_id);
        }
        auto ofastq = TextWriter::create(fastq_output_path, output_compression);
        auto omsa = TextWriter::create(msa_output_path, output_compression);
        std::array<char, 1 << 15> copy_buffer{};
        auto append_file_to_writer = [&](const std::string& path, TextWriter& writer) {
            std::ifstream input(path, std::ios::binary);
            if (!input.is_open()) {
                std::cout << "Failed to open temporary file " << path << "\n";
                exit(-1);
            }
            while (input) {
                input.read(copy_buffer.data(), static_cast<std::streamsize>(copy_buffer.size()));
                std::streamsize count = input.gcount();
                if (count > 0) {
                    writer.write(copy_buffer.data(), static_cast<size_t>(count));
                }
            }
        };
        for (size_t thread_id = 0; thread_id < thread_count; thread_id++) {
            threads[thread_id].join();

            std::string fastq_t_filename = o_filename_prefix + "fastq" + std::to_string(thread_id);
            append_file_to_writer(fastq_t_filename, *ofastq);
            remove(fastq_t_filename.c_str());

            std::string msa_t_filename = o_filename_prefix + "msa" + std::to_string(thread_id);
            append_file_to_writer(msa_t_filename, *omsa);
            remove(msa_t_filename.c_str());
        }
    }
}

int main(int argc, char** argv) {
    parse_flags(argc, argv);
    run_consensus();
    return 0;
}
