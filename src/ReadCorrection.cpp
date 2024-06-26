// ReadCorrection.cpp
#include "ReadCorrection.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <iostream>
#include <vector>

using namespace boost;

ReadCorrection::ReadCorrection(Graph graph, read_arguments args, std::unordered_map<std::vector<seqan3::dna5>, Vertex, std::hash<std::vector<seqan3::dna5>>> read2vertex, std::unordered_map<Vertex, std::vector<seqan3::dna5>, std::hash<Vertex>> vertex2read) : graph_(std::move(graph)), args(args), read2vertex_(std::move(read2vertex)), vertex2read_(std::move(vertex2read)) {}

bool ReadCorrection::is_isolated(Vertex v) {
    return boost::degree(v, graph_) == 0;
}

std::vector<Vertex> ReadCorrection::get_connected_nodes_with_weight(Vertex v, int w) {
    std::vector<Vertex> connected_nodes;
    for (auto edge : boost::make_iterator_range(boost::out_edges(v, graph_))) {
        if (graph_[edge].weight == w) {
            Vertex target = boost::target(edge, graph_);
            connected_nodes.push_back(target);
        }
    }
    return connected_nodes;
}

std::vector<std::vector<seqan3::dna5>> ReadCorrection::get_low_frequency_reads(const std::map<std::vector<seqan3::dna5>, uint32_t>& read2count) {
    std::vector<std::vector<seqan3::dna5>> low_frequency_reads;
    for (const auto& pair : read2count) {
        if (pair.second < args.freq_thresh) {
            low_frequency_reads.push_back(pair.first);
        }
    }
    return low_frequency_reads;
}

void ReadCorrection::correction_isolates(const std::map<std::vector<seqan3::dna5>, uint32_t>& read2count){
    std::vector<Vertex> isolates;
    for (const auto& pair : read2count) {
        auto cur_vertex = read2vertex_[pair.first];
        if (is_isolated(cur_vertex)){
            isolates.push_back(cur_vertex);
        }
    }

}

void ReadCorrection::correction_process(const std::vector<std::vector<seqan3::dna5>>& low_count_reads, int w) {
    for (const auto& read : low_count_reads) {
        auto cur_vertex = read2vertex_[read];

        if (!is_isolated(cur_vertex)) {
            auto connected_nodes = get_connected_nodes_with_weight(cur_vertex, w);
            int edge_num = connected_nodes.size();

            if (edge_num == 1) {
                Vertex neighbor = connected_nodes[0];
                if (graph_[neighbor].count > args.freq_thresh) {
                    graph_[neighbor].count += graph_[cur_vertex].count;
                    graph_[neighbor].add_ids(graph_[cur_vertex].ids);
                    boost::clear_vertex(cur_vertex, graph_);
                    boost::remove_vertex(cur_vertex, graph_);
                }
            } else if (edge_num >= 2) {
                std::vector<Vertex> eligible_nodes;
                int total_count = 0;
                for (Vertex node : connected_nodes) {
                    if (graph_[node].count > args.freq_thresh) {
                        eligible_nodes.push_back(node);
                        total_count += graph_[node].count;
                    }
                }

                if (eligible_nodes.size() == 1) {
                    Vertex target = eligible_nodes[0];
                    graph_[target].count += graph_[cur_vertex].count;
                    graph_[target].add_ids(graph_[cur_vertex].ids);
                    boost::clear_vertex(cur_vertex, graph_);
                    boost::remove_vertex(cur_vertex, graph_);
                } else if (!eligible_nodes.empty()) {
                    // Sort eligible_nodes by their count in descending order
                    std::sort(eligible_nodes.begin(), eligible_nodes.end(), 
                              [this](Vertex a, Vertex b) {
                                  return graph_[a].count > graph_[b].count;
                              });

                    int original_count = graph_[cur_vertex].count;
                    int remaining_count = original_count;

                    for (Vertex node : eligible_nodes) {
                        if (remaining_count == 0) break;

                        double proportion = static_cast<double>(graph_[node].count) / total_count;
                        int added_count = std::round(remaining_count * proportion);

                        graph_[node].count += added_count;
                        // add ids to the node
                        std::vector<std::string> first_m_ids(graph_[cur_vertex].ids.begin(), graph_[cur_vertex].ids.begin() + added_count);
                        graph_[cur_vertex].ids.erase(graph_[cur_vertex].ids.begin(), graph_[cur_vertex].ids.begin() + added_count);
                        graph_[node].add_ids(first_m_ids);

                        remaining_count -= added_count;
                    }

                    // Assign any remaining count to the node with the highest count
                    if (remaining_count > 0) {
                        Vertex highest_count_node = eligible_nodes.front(); // First node has highest count
                        
                        graph_[highest_count_node].count += remaining_count;
                        graph_[highest_count_node].add_ids(graph_[cur_vertex].ids);
                    }

                    boost::clear_vertex(cur_vertex, graph_);
                    boost::remove_vertex(cur_vertex, graph_);
                }
            }
        }
    }
}

void ReadCorrection::correction_main(const std::map<std::vector<seqan3::dna5>, uint32_t>& read2count){
    auto low_count_reads = get_low_frequency_reads(read2count);
    for (int w_i = 1; w_i <= args.max_edit_dis; w_i++){
        correction_process(low_count_reads, w_i);
    }
    // correction_isolates(read2count);
}

// std::string ReadCorrection::get_output_filename(std::string output_type) {
//     std::filesystem::path input_path(args.input_data);
//     std::filesystem::path output_path(args.output_dir);
//     std::string output_filename = input_path.stem().string() + output_type + input_path.extension().string();
//     return (output_path / output_filename).string();
// }

// bool ReadCorrection::is_fastq_by_extension(const std::string &input_file) {
//     std::string ext = std::filesystem::path(input_file).extension().string();
//     std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower); // Convert to lower case for case-insensitive comparison

//     if (ext == ".gz") {
//         std::string stem_ext = std::filesystem::path(input_file).stem().extension().string();
//         std::transform(stem_ext.begin(), stem_ext.end(), stem_ext.begin(), ::tolower); // Convert to lower case for case-insensitive comparison
        
//         if (stem_ext == ".fastq" || stem_ext == ".fq") {
//             return true;
//         } else if (stem_ext == ".fasta" || stem_ext == ".fa") {
//             return false;
//         } else {
//             throw std::runtime_error("Unknown file extension before .gz: " + stem_ext);
//         }
//     } else {
//         if (ext == ".fastq" || ext == ".fq") {
//             return true;
//         } else if (ext == ".fasta" || ext == ".fa") {
//             return false;
//         } else {
//             throw std::runtime_error("Unknown file extension: " + ext);
//         }
//     }
// }

// void ReadCorrection::write_graph2dataset(std::map<std::string, std::vector<seqan3::phred42>>, ){
//     // Define the type of sequence file input
//     seqan3::sequence_file_input input{args.input_data};

//     auto output_correction_file = get_output_filename(std::string ".correction");
//     std::ofstream output_stream(output_correction_file);

//     auto is_fastq = is_fastq_by_extension(args.input_data);

//     // Define the type of sequence file output
//     if (is_fastq) {
//         seqan3::sequence_file_output output{output_stream, seqan3::format_fastq{}};
//         for (const auto &vertex : boost::make_iterator_range(vertices(graph_))) {
//             const auto &props = graph_[vertex];
//             for (const auto &id : props.ids) {
//                 output.emplace_back(props.read, id + " " + id2description_[id], id2quality_[id]);
//             }
//         }
//     } else {
//         seqan3::sequence_file_output output{output_stream, seqan3::format_fasta{}};
//         for (const auto &vertex : boost::make_iterator_range(vertices(graph_))) {
//             const auto &props = graph_[vertex];
//             for (const auto &id : props.ids) {
//                 output.emplace_back(props.read, id + " " + id2description_[id]);
//             }
//         } 
//     }
// }
