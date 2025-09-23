//
// Created by borabi on 20/12/17.
//

#include "cluster.h"

// #include <pthread.h>
#include <thread>
#include <mutex>
#include <stack>
#include <algorithm>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <exception>
#include <iterator>
#include <utility>
// Debug includes
#include <sstream>
#include <iomanip>

#include "text_io.h"

using namespace std;


// extern variables declarations
cluster_id_t cluster_count = 0;
node_id_to_cluster_id_vector node_to_cluster_vector;

// locally global variables
#define ASCII_SIZE 256
#define MAX_TMP_FILE_COUNT 25.0
bool valid_base [ASCII_SIZE];
node_id_to_node_id_vector_of_vectors* graph_ptr;
mutex graph_lock;
vector<vector<bool> > all_masks;
long max_memory_use = 0;


void cluster(){
    time_t start;
    graph_ptr = new node_id_to_node_id_vector_of_vectors();
    if (!silent) {
        cout << "Adding edges due to barcode barcode similarity\n";
    }
    start = time(NULL);
    barcode_similarity();
    if (!silent) {
        cout << "Adding edges due to barcodes similarity took: " << difftime(time(NULL), start) << "\n";
    }
    if (!silent && print_mem){
        cout << "Memory after adding edges:\n\t" << get_memory_use() << "MB\n";
    }

    if (!silent) {
        cout << "Extracting clusters\n";
    }
    start = time(NULL);
    extract_clusters();
    if (!silent) {
        cout << "Extracting clusters took: " << difftime(time(NULL), start) << "\n";
    }
    if (!silent && print_mem){
        cout << "Memory extracting clusters:\n\t" << get_memory_use() << "MB\n";
    }
    delete graph_ptr;
    if (!silent && print_mem){
        cout << "Memory after releasing graph:\n\t" << get_memory_use() << "MB\n";
    }

    if (!silent) {
        cout << "Outputting clusters\n";
    }
    start = time(NULL);
    output_clusters();
    if (!silent) {
        cout << "Outputting clusters took: " << difftime(time(NULL), start) << "\n";
    }

}

void process_lsh(masked_barcode_to_barcode_id_unordered_map* lsh_ptr,
                    node_id_to_node_id_vector_of_vectors* local_graph_ptr) {
    for (auto kv = (*lsh_ptr).begin(); kv!= (*lsh_ptr).end(); kv++) {
        for (barcode_id_t bid: kv->second){
            for (node_id_t node : (*barcode_to_nodes_vector_ptr)[bid]) {
                for (barcode_id_t bid_o: kv->second){
                    if (bid == bid_o){
                        continue;
                    }
                    vector<node_id_t> good_neighbors = get_good_neighbors(node, (*barcode_to_nodes_vector_ptr)[bid_o]);
                    vector<node_id_t> result;
                    set_union((*local_graph_ptr)[node].begin(), (*local_graph_ptr)[node].end(),
                        good_neighbors.begin(), good_neighbors.end(),
                        back_inserter(result)
                    );
                    (*local_graph_ptr)[node] = std::move(result);
                }
            }
        }
    }
}

void process_identical_barcode_nodes(uint8_t barcode_id_remainder) {
    if (!silent) {
        stringstream stream;
        stream << "Adding edges between nodes of identical barcodes with thread " << (int)barcode_id_remainder << "\n";
        cout << stream.str();
    }
    for (barcode_id_t i = barcode_id_remainder; i < barcode_count; i+=thread_count) {
        for (node_id_t node : (*barcode_to_nodes_vector_ptr)[i]) {
            vector<node_id_t> good_neighbors = get_good_neighbors(node, (*barcode_to_nodes_vector_ptr)[i]);
            vector<node_id_t> result;
            set_union((*graph_ptr)[node].begin(), (*graph_ptr)[node].end(),
                      good_neighbors.begin(), good_neighbors.end(),
                      back_inserter(result)
                      );
        (*graph_ptr)[node] = std::move(result);
        }
    }
}

void lsh_mask(size_t mask_remainder) {
    time_t start;
    time_t build_time = 0, process_time = 0;
    node_id_to_node_id_vector_of_vectors local_graph(node_count);
    char masked_barcode_buffer[150];
    masked_barcode_buffer[barcode_length_1+barcode_length_2-error_tolerance] = '\0';
    stringstream stream;
    for (size_t i = 0; i < all_masks.size(); i++) {
        if (i % thread_count != mask_remainder) {
            continue;
        }
        start = time(NULL);
        vector<bool> mask = all_masks[i];
        masked_barcode_to_barcode_id_unordered_map lsh;
        if (!silent) {
            string current_mask_bin;
            for (bool p: mask) {
                current_mask_bin += p ? "1" : "0";
            }
            stream = stringstream();
            stream << current_mask_bin << " is assigned to thread "<< mask_remainder << "\n";
            cout << stream.str();
        }
        string masked_barcode;
        for (barcode_id_t i = 0; i < barcode_count; i++) {
            masked_barcode = mask_barcode((*barcodes_ptr)[i], mask, masked_barcode_buffer);
            if (masked_barcode != "0"){
                lsh[masked_barcode].push_back(i);
            }
        }
        build_time += difftime(time(NULL), start);
        if (!silent) {
            stream = stringstream();
            stream << "Thread "<< mask_remainder <<" built LSH in: " << difftime(time(NULL), start) << "\n";
            cout << stream.str();
        }
        start = time(NULL);
        process_lsh(&lsh, &local_graph);
        process_time += difftime(time(NULL), start);
        if (!silent) {
            stream = stringstream();
            stream << "Thread "<< mask_remainder <<" processed LSH in: " << difftime(time(NULL), start) << "\n";
            cout << stream.str();
        }
    }

    if (!silent) {
        stream = stringstream();
        stream << "On thread " << mask_remainder << " building all LSH took: " << build_time << "\n";
        stream << "On thread " << mask_remainder << " processing all LSH took: " << process_time << "\n";
        cout << stream.str();
    }
    start = time(NULL);
    if (!silent) {
        stream = stringstream();
        stream << "On thread " << mask_remainder << " merging local graph with global graph\n";
        cout << stream.str();
    }
    merge_graphs(&local_graph);
    node_id_to_node_id_vector_of_vectors().swap(local_graph);
    if (!silent) {
        stream = stringstream();
        stream << "On thread " << mask_remainder << " merging took " << difftime(time(NULL), start) <<"\n";
        cout << stream.str();
    }

}

void merge_graphs(node_id_to_node_id_vector_of_vectors* local_graph_ptr) {
    graph_lock.lock();
    if (sort_clusters) {
        max_memory_use = max((long) (get_memory_use()*1024), (long) max_memory_use);
    }
    if ((*graph_ptr).size() == 0) {
        (*graph_ptr) = std::move(*local_graph_ptr);
        graph_lock.unlock();
        return;
    }
    for (node_id_t node = 0; node < node_count; node++) {
        vector<node_id_t> result;
        set_union((*local_graph_ptr)[node].begin(), (*local_graph_ptr)[node].end(),
                    (*graph_ptr)[node].begin(), (*graph_ptr)[node].end(),
                    back_inserter(result)
        );
        (*graph_ptr)[node] = std::move(result);
    }
    graph_lock.unlock();
}

void barcode_similarity(){
    std::fill(std::begin(valid_base), std::end(valid_base), false);
    auto mark_base = [&](char base) {
        valid_base[static_cast<unsigned char>(base)] = true;
    };
    mark_base('A');
    mark_base('C');
    mark_base('G');
    mark_base('T');
    mark_base('a');
    mark_base('c');
    mark_base('g');
    mark_base('t');

    size_t mask_count = 1;
    vector<bool> mask(barcode_length_1+barcode_length_2, false);
    std::fill(mask.begin() + error_tolerance, mask.end(), true);
    for (int i = barcode_length_1+barcode_length_2; i > barcode_length_1+barcode_length_2 - error_tolerance; i--) {
        mask_count *= i;
        mask_count /= barcode_length_1+barcode_length_2 - i + 1;
    }
    all_masks.reserve(mask_count);
    do{
        all_masks.push_back(mask);

    } while (std::next_permutation(mask.begin(), mask.end()));
    if (!silent) {
        cout << "Number of masks is " << all_masks.size() << "\n";
    }

    thread* thread_array = NULL;
    time_t start = time(NULL);
    if (thread_count > 1) {
        thread_array = new thread[thread_count];
        for (int t_id = 0; t_id < thread_count; t_id++) {
            thread_array[t_id] = thread(lsh_mask, t_id);
            stringstream stream;
            stream << "Created thread " << t_id << "\n";
            cout << stream.str();
        }
        for (int t_id = 0; t_id < thread_count; t_id++) {
            thread_array[t_id].join();
            stringstream stream;
            stream << "Joined thread " << t_id << "\n";
            cout << stream.str();
        }
    } else {
        lsh_mask(0);
    }
    if (!silent) {
        cout << "Building the graph on " << thread_count << " thread(s) took " <<  difftime(time(NULL), start) << "\n";
    }
    // barcodes are no longer needed
    delete barcodes_ptr;
    if (thread_count > 1) {
        for (int t_id = 0; t_id < thread_count; t_id++) {
            thread_array[t_id] = thread(process_identical_barcode_nodes, t_id);
        }
        for (int t_id = 0; t_id < thread_count; t_id++) {
            thread_array[t_id].join();
            stringstream stream;
            stream << "Joined thread " << t_id << "\n";
            cout << stream.str();
        }
        delete [] thread_array;
    } else {
        process_identical_barcode_nodes(0);
    }
    // barcode id to node id's is no longer needed
    delete barcode_to_nodes_vector_ptr;
    delete node_to_minimizers_ptr;
}

string mask_barcode(const string& barcode, const vector<bool>& mask, char* masked_barcode_buffer){
    int pos = 0;
    for (int i = 0; i < barcode_length_1+barcode_length_2; i++) {
        if (mask[i]) {
            if (valid_base[(uint8_t) barcode.at(i)] == false) {
                return "0";
            }
            masked_barcode_buffer[pos] = barcode.at(i);
            pos++;
        }
    }
    return string(masked_barcode_buffer);
}

vector<node_id_t> get_good_neighbors(node_id_t node, const vector<node_id_t> &neighbors){
    vector<node_id_t> good_neighbors;
    for (node_id_t neighbor : neighbors) {
        if (node != neighbor && !unmatched_minimimizers(node, neighbor)) {
            good_neighbors.push_back(neighbor);
        }
    }
    return good_neighbors;
}

bool unmatched_minimimizers(node_id_t node_id, node_id_t neighbor_id){
    int matched_minimimizers_1 = 0;
    int matched_minimimizers_2 = 0;
    for (int i =0; i < minimizer_count; i++) {
        matched_minimimizers_1 += (*node_to_minimizers_ptr)[node_id][i] == (*node_to_minimizers_ptr)[neighbor_id][i];
        matched_minimimizers_2 += (*node_to_minimizers_ptr)[node_id][i+minimizer_count] == (*node_to_minimizers_ptr)[neighbor_id][i+minimizer_count];
    }
    return !(matched_minimimizers_1 >= minimizer_threshold && matched_minimimizers_2 >= minimizer_threshold);
}

void extract_clusters(){
    vector<bool> pushed(node_count, false);
    stack<node_id_t> opened;
    node_to_cluster_vector.reserve(node_count);
    cluster_count = 0;

    for (node_id_t node = 0; node < node_count; node++) {
        if (!pushed[node]) {
            opened.push(node);
            pushed[node] = true;
            while(!opened.empty()) {
                node_id_t current_node = opened.top();
                opened.pop();
                node_to_cluster_vector[current_node] = cluster_count;
                for (node_id_t neighbor: (*graph_ptr)[current_node]) {
                    if (!pushed[neighbor]) {
                        opened.push(neighbor);
                        pushed[neighbor] = true;
                        node_to_cluster_vector[neighbor] = cluster_count;
                    }
                }
            }
            cluster_count++;
        }
    }
}

void output_clusters(){
    FastqRecord record_1;
    FastqRecord record_2;
    try {
        FastqReader fastq_reader_1(input_1, input_compression);
        FastqReader fastq_reader_2(input_2, input_compression);

        auto build_cluster_record = [&](node_id_t current_read_node,
                                        read_id_t current_read,
                                        const std::string& name_1,
                                        const std::string& sequence_1,
                                        const std::string& quality_1,
                                        const std::string& name_2,
                                        const std::string& sequence_2,
                                        const std::string& quality_2) {
            std::string record;
            if (compact_cluster_output) {
                record.reserve(32);
                record += std::to_string(node_to_cluster_vector[current_read_node]);
                record.push_back('\t');
                record += std::to_string(current_read_node);
                record.push_back('\t');
                record += std::to_string(current_read);
            } else {
                record.reserve(name_1.size() + name_2.size() + sequence_1.size() + sequence_2.size() + quality_1.size() + quality_2.size() + 64);
                record += std::to_string(node_to_cluster_vector[current_read_node]);
                record.push_back('\t');
                record += std::to_string(current_read_node);
                record.push_back('\t');
                record += std::to_string(current_read);
                record.push_back('\t');
                record += name_1;
                record.push_back('\t');
                record += sequence_1;
                record.push_back('\t');
                record += quality_1;
                record.push_back('\t');
                record += name_2;
                record.push_back('\t');
                record += sequence_2;
                record.push_back('\t');
                record += quality_2;
            }
            record.push_back('\n');
            return record;
        };

        const std::string cluster_extension = compressionTypeExtension(cluster_output_compression);
        const std::string cluster_output_path = output_prefix + "cluster" + cluster_extension;

        if (!sort_clusters) {
            read_id_t current_read = 0;
            auto cluster_writer = TextWriter::create(cluster_output_path, cluster_output_compression);

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

                if (record_1.sequence.size() != record_1.quality.size() ||
                    record_2.sequence.size() != record_2.quality.size()) {
                    cerr << "ERROR: Sequence and quality length mismatch detected.\n";
                    cerr << "name_1\t" << record_1.name << "\n";
                    cerr << "sequence_1\t" << record_1.sequence << "\n";
                    cerr << "quality_1\t" << record_1.quality << "\n";
                    cerr << "name_2\t" << record_2.name << "\n";
                    cerr << "sequence_2\t" << record_2.sequence << "\n";
                    cerr << "quality_2\t" << record_2.quality << "\n";
                    exit(-1);
                }

                node_id_t current_read_node = read_to_node_vector[current_read];
                std::string record = build_cluster_record(current_read_node,
                                                          current_read,
                                                          record_1.name,
                                                          record_1.sequence,
                                                          record_1.quality,
                                                          record_2.name,
                                                          record_2.sequence,
                                                          record_2.quality);
                cluster_writer->write(record);
                current_read++;
        }
        return;
    }

        read_id_t current_read = 0;
        std::unique_ptr<TextWriter> cluster_writer;
        if (max_memory_use == 0) {
            max_memory_use = 1024;
        }
        size_t min_records_per_tmp_file = max_memory_use/4;
        cout << "min_records_per_tmp_file " << min_records_per_tmp_file << "\n";
        size_t temp_out_count;
        cout << "There are " << cluster_count << " clusters\n";
        temp_out_count = (unsigned long) ceil(float(read_count)/float(min_records_per_tmp_file));
        cout << "There are " << temp_out_count << " temp files\n";
        temp_out_count = min(MAX_TMP_FILE_COUNT, (double) temp_out_count);
        cout << "There are " << temp_out_count << " temp files\n";
        vector<ofstream> temp_out_files(temp_out_count);
        vector<string> temp_out_names(temp_out_count);
        for (size_t i = 0; i < temp_out_count; i++) {
            temp_out_names[i] = output_prefix + "temp_" + to_string(i);
            temp_out_files[i] = ofstream(temp_out_names[i]);
        }

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

            if (record_1.sequence.size() != record_1.quality.size() ||
                record_2.sequence.size() != record_2.quality.size()) {
                cerr << "ERROR: Sequence and quality length mismatch detected.\n";
                cerr << "name_1\t" << record_1.name << "\n";
                cerr << "sequence_1\t" << record_1.sequence << "\n";
                cerr << "quality_1\t" << record_1.quality << "\n";
                cerr << "name_2\t" << record_2.name << "\n";
                cerr << "sequence_2\t" << record_2.sequence << "\n";
                cerr << "quality_2\t" << record_2.quality << "\n";
                exit(-1);
            }

            node_id_t current_read_node = read_to_node_vector[current_read];
            size_t current_temp_out_id = node_to_cluster_vector[current_read_node] % temp_out_count;
            std::string record = build_cluster_record(current_read_node,
                                                      current_read,
                                                      record_1.name,
                                                      record_1.sequence,
                                                      record_1.quality,
                                                      record_2.name,
                                                      record_2.sequence,
                                                      record_2.quality);
            temp_out_files[current_temp_out_id] << record;
            current_read++;
        }
        read_id_to_node_id_vector().swap(read_to_node_vector);
        node_id_to_cluster_id_vector().swap(node_to_cluster_vector);

        cluster_writer = TextWriter::create(cluster_output_path, cluster_output_compression);
        for (size_t i = 0; i < temp_out_count; i++) {
            cout << "Processing file " << temp_out_names[i] << "\n";
            temp_out_files[i].close();
            ifstream temp_file;
            temp_file.open(temp_out_names[i]);
            vector<string> records(size_t(ceil((double)cluster_count/(double)temp_out_count)), "");
            string record;
            while(getline(temp_file, record)) {
                size_t cluster_id = stoi(record.substr(0, record.find("\t"))) % temp_out_count;
                records[cluster_id]+= record + "\n";
            }
            for (const string& record_line : records) {
                if (!record_line.empty()) {
                    cluster_writer->write(record_line);
                }
            }
            temp_file.close();
            remove(temp_out_names[i].c_str());
        }
    } catch (const std::exception& ex) {
        cout << ex.what() << "\n";
        exit(-1);
    }
}
