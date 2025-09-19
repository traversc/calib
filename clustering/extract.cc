//
// Created by borabi on 19/12/17.
//

#include "extract.h"
#include "global.h"
#include <algorithm>
#include <exception>
#include <utility>

// Debug includes
#include <sstream>
#include <iomanip>

using namespace std;

#define ASCII_SIZE 128
#define BYTE_SIZE 8
#define ENCODE_SIZE 2

node_id_t node_count = 0;
read_id_t read_count = 0;
barcode_id_t barcode_count = 0;

// node_to_read_id finds reads with identical barcodes and minimizers
typedef unordered_map<NodePtr, vector<read_id_t>, NodePtrHash, NodePtrEqual> node_to_read_id_unordered_map;
// barcode_to_node_id finds nodes with identical barcodes
typedef unordered_map<barcode_t, vector<node_id_t> > barcode_to_node_id_unordered_map;

read_id_to_node_id_vector read_to_node_vector;
// // These data structures are accessed by different threads; create them on heap
barcode_vector* barcodes_ptr;
barcode_id_to_node_ids_vector* barcode_to_nodes_vector_ptr;
node_id_to_minimizers_vector* node_to_minimizers_ptr;

// Varialble to encode DNA nucleotides into 2-bits and keeping track
// of segments that don't have valid kmers (all kmers have >= 1 N)
minimizer_t encode [ASCII_SIZE];
minimizer_t invalid_kmer = (minimizer_t) -1;

string minimizer_t_to_dna(minimizer_t minimizer, size_t size) {
    string s = "";
    for (int i = 0; i < (int) size; i++) {
        switch (minimizer & 0x3) {
            case 0x0:
                s = "A" + s;
                break;
            case 0x1:
                s = "C" + s;
                break;
            case 0x2:
                s = "G" + s;
                break;
            case 0x3:
                s = "T" + s;
                break;
        }
        minimizer = minimizer >> ENCODE_SIZE;
    }
    return s;
}

void extract_barcodes_and_minimizers() {
    // Create those on the heap
    barcodes_ptr = new barcode_vector();
    barcode_to_nodes_vector_ptr = new barcode_id_to_node_ids_vector();
    node_to_minimizers_ptr = new node_id_to_minimizers_vector();

    node_to_read_id_unordered_map node_to_read_map;

    for (int i = 0; i < ASCII_SIZE; i++) {
        encode[i] = (minimizer_t) -1;
    }
    auto set_encode = [&](char base, minimizer_t value) {
        encode[static_cast<unsigned char>(base)] = value;
    };
    set_encode('A', 0x0);
    set_encode('C', 0x1);
    set_encode('G', 0x2);
    set_encode('T', 0x3);
    set_encode('a', 0x0);
    set_encode('c', 0x1);
    set_encode('g', 0x2);
    set_encode('t', 0x3);

    read_count = 0;
    if (!silent && print_mem){
        cout << "Memory before reading FASTQ:\n\t" << get_memory_use() << "MB\n";
    }
    // Processing FASTQ files one read at a time
    FastqRecord record_1;
    FastqRecord record_2;
    try {
        FastqReader fastq_reader_1(input_1, input_compression);
        FastqReader fastq_reader_2(input_2, input_compression);

        while (true) {
            bool has_first = fastq_reader_1.readRecord(record_1);
            bool has_second = fastq_reader_2.readRecord(record_2);
            if (!has_first || !has_second) {
                if (has_first != has_second) {
                    cout << "Input FASTQ files have different numbers of reads.\n";
                    exit(-1);
                }
                break;
            }

            string& sequence_1 = record_1.sequence;
            string& sequence_2 = record_2.sequence;

            int s1_length = static_cast<int>(sequence_1.size());
            int s2_length = static_cast<int>(sequence_2.size());

            NodePtr current_node_ptr = new Node();
            // Extracting the barcode from the start of both mates
            if (s1_length >= barcode_length_1 && s2_length >= barcode_length_2) {
                current_node_ptr->barcode =
                    sequence_1.substr(0, barcode_length_1) +
                    sequence_2.substr(0, barcode_length_2);
            } else {
                current_node_ptr->barcode = string (barcode_length_1 + barcode_length_2, 'N');
            }

            s1_length -= barcode_length_1;
            s2_length -= barcode_length_2;

            s1_length -= ignored_sequence_prefix_length;
            s2_length -= ignored_sequence_prefix_length;

            // Splitting the remaining sequence into ~ equally sized segments, and extracting minimizers from each
            int s1_seg_length = minimizer_count > 0 ? s1_length / minimizer_count : 0;
            if (s1_seg_length >= kmer_size) {
                int start = barcode_length_1 + ignored_sequence_prefix_length;
                for (int i = 0; i < minimizer_count; i++) {
                    current_node_ptr->minimizers[i] = minimizer(sequence_1, start, s1_seg_length);
                    start += s1_seg_length;
                }
            } else {
                for (int i = 0; i < minimizer_count; i++) {
                    current_node_ptr->minimizers[i] = --invalid_kmer;
                }
            }

            int s2_seg_length = minimizer_count > 0 ? s2_length / minimizer_count : 0;
            if (s2_seg_length >= kmer_size) {
                int start = barcode_length_2 + ignored_sequence_prefix_length;
                for (int i = minimizer_count; i < minimizer_count*2; i++) {
                    current_node_ptr->minimizers[i] = minimizer(sequence_2, start, s2_seg_length);
                    start += s2_seg_length;
                }
            } else {
                for (int i = minimizer_count; i < minimizer_count*2; i++) {
                    current_node_ptr->minimizers[i] = --invalid_kmer;
                }
            }

            if (node_to_read_map.find(current_node_ptr) != node_to_read_map.end()) {
                node_to_read_map[current_node_ptr].push_back(read_count);
                delete current_node_ptr;
            } else {
                node_to_read_map.emplace(current_node_ptr, vector<read_id_t>{read_count});
            }
            read_count++;
        }
    } catch (const std::exception& ex) {
        cout << ex.what() << "\n";
        exit(-1);
    }
    if (!silent && print_mem){
        cout << "Memory right after reading FASTQ:\n\t" << get_memory_use() << "MB\n";
    }

    read_to_node_vector.reserve(read_count);
    (*node_to_minimizers_ptr).reserve(node_to_read_map.size());

    if (!silent && print_mem){
        cout << "Memory after reserving for read_to_node_vector & node_to_minimizers:\n\t" << get_memory_use() << "MB\n";
    }

    node_count = 0;
    barcode_to_node_id_unordered_map barcode_to_node_map;
    for (auto kv : node_to_read_map) {
        NodePtr current_node_ptr = kv.first;
        (*node_to_minimizers_ptr).emplace_back(std::move(current_node_ptr->minimizers));
        barcode_to_node_map[current_node_ptr->barcode].push_back(node_count);
        delete current_node_ptr;
        for (read_id_t rid : kv.second) {
            read_to_node_vector[rid] = node_count;
        }
        node_count++;
    }
    if (!silent && print_mem){
        cout << "Memory after filling barcode_to_node_map:\n\t" << get_memory_use() << "MB\n";
    }
    node_to_read_id_unordered_map().swap(node_to_read_map);
    if (!silent && print_mem){
        cout << "Memory after releasing node_to_read_map:\n\t" << get_memory_use() << "MB\n";
    }

    barcode_count = barcode_to_node_map.size();
    (*barcodes_ptr).reserve(barcode_count);
    (*barcode_to_nodes_vector_ptr).reserve(barcode_count);
    if (!silent && print_mem){
        cout << "Memory after reserving barcode_to_nodes_vector:\n\t" << get_memory_use() << "MB\n";
    }

    for (auto kv : barcode_to_node_map) {
        (*barcodes_ptr).emplace_back(std::move(kv.first));
        sort(kv.second.begin(), kv.second.end());
        (*barcode_to_nodes_vector_ptr).emplace_back(std::move(kv.second));
    }
    if (!silent && print_mem){
        cout << "Memory after filling barcodes & barcode_to_nodes_vector:\n\t" << get_memory_use() << "MB\n";
    }
    barcode_to_node_id_unordered_map().swap(barcode_to_node_map);
    if (!silent && print_mem){
        cout << "Memory after releasing barcode_to_node_map:\n\t" << get_memory_use() << "MB\n";
    }

    if (!silent) {
        cout << "Read count: " << read_count << "\n";
        cout << "Node count: " << node_count << "\n";
        cout << "Barcode count: " << barcode_count << "\n";
    }
}

// Extract the lexicographically minimum k-mer in a given string with start and range
minimizer_t minimizer(string& seq, int start, int length){
    int end = start + length;
    minimizer_t current_k_mer = (minimizer_t) -1;
    // Biggest possible k-mer is all 1's
    minimizer_t min_k_mer = (minimizer_t) -1;
    // Ignoring leftmost bits depending on k-mer size
    minimizer_t kmer_size_mask = (minimizer_t) -1;
    // TODO: Compute this once?
    kmer_size_mask >>= (sizeof(minimizer_t)*BYTE_SIZE-kmer_size*ENCODE_SIZE);

    if (start + kmer_size > end) {
        return --invalid_kmer;
    }

    // Building the first k-mer
    for (int i = start; i < start + kmer_size && i < end; i++) {
        current_k_mer <<= ENCODE_SIZE;
        current_k_mer |= encode[(size_t)seq.at(i)];
        // Hit a non ACTG nucleotide
        if (encode[(size_t)seq.at(i)] == (minimizer_t) -1) {
            current_k_mer = (minimizer_t) -1;
            start = i+1;
            // Hit yet another non nucleotide
            if (start + kmer_size > end) {
                return --invalid_kmer;
            }
        }
    }
    current_k_mer &= kmer_size_mask;
    min_k_mer = (min_k_mer < current_k_mer)  ? min_k_mer : current_k_mer;

    for (int i = start + kmer_size; i < end; i++) {
        current_k_mer <<= ENCODE_SIZE;
        current_k_mer &= kmer_size_mask;
        current_k_mer |= encode[(size_t)seq.at(i)];
        // Re-building the first k-mer if we hit a non ACTG nucleotide
        if (encode[(size_t)seq.at(i)] == (minimizer_t) -1) {
            current_k_mer = (minimizer_t) -1;
            i++;
            int j;
            for (j = i; j < i + kmer_size && j < end; j++) {
                current_k_mer <<= ENCODE_SIZE;
                current_k_mer |= encode[(size_t)seq.at(j)];

                // Hit a non ACTG nucleotide
                if (encode[(size_t)seq.at(j)] == (minimizer_t) -1) {
                    i = j+1;
                }
            }
            current_k_mer &= kmer_size_mask;
            i=j;
        }
        min_k_mer = (min_k_mer < current_k_mer) ? min_k_mer : current_k_mer;
    }
    return min_k_mer != (minimizer_t) -1 ? min_k_mer : --invalid_kmer;
}
