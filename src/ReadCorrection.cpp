// ReadCorrection.cpp
#include "ReadCorrection.hpp"

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

void ReadCorrection::correction_process(const std::vector<std::vector<seqan3::dna5>>& low_count_reads, int w) {
    std::vector<Vertex> isolated_node_ids;

    for (const auto& read : low_count_reads) {
        auto cur_id = read2vertex_[read];

        if (is_isolated(cur_id)) {
            isolated_node_ids.push_back(cur_id);
        } else {
            auto connected_nodes = get_connected_nodes_with_weight(cur_id, w);
            int edge_num = connected_nodes.size();

            if (edge_num == 1) {
                Vertex neighbor = connected_nodes[0];
                if (graph_[neighbor].count > args.freq_thresh) {
                    graph_[neighbor].count += graph_[cur_id].count;
                    boost::clear_vertex(cur_id, graph_);
                    boost::remove_vertex(cur_id, graph_);
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
                    graph_[target].count += graph_[cur_id].count;
                    boost::clear_vertex(cur_id, graph_);
                    boost::remove_vertex(cur_id, graph_);
                } else if (!eligible_nodes.empty()) {
                    // Sort eligible_nodes by their count in descending order
                    std::sort(eligible_nodes.begin(), eligible_nodes.end(), 
                              [this](Vertex a, Vertex b) {
                                  return graph_[a].count > graph_[b].count;
                              });

                    int original_count = graph_[cur_id].count;
                    int remaining_count = original_count;

                    for (Vertex node : eligible_nodes) {
                        if (remaining_count == 0) break;

                        double proportion = static_cast<double>(graph_[node].count) / total_count;
                        int added_count = std::round(remaining_count * proportion);

                        graph_[node].count += added_count;
                        remaining_count -= added_count;
                    }

                    // Assign any remaining count to the node with the highest count
                    if (remaining_count > 0) {
                        Vertex highest_count_node = eligible_nodes.front(); // First node has highest count
                        graph_[highest_count_node].count += remaining_count;
                    }

                    boost::clear_vertex(cur_id, graph_);
                    boost::remove_vertex(cur_id, graph_);
                }
            }
        }
    }

    // Handle isolated nodes if needed
    // For example, you might want to log or process isolated_node_ids further.
}

void ReadCorrection::correction_main(const std::map<std::vector<seqan3::dna5>, uint32_t>& read2count){
    auto low_count_reads = get_low_frequency_reads(read2count);
    for (int w_i = 1; w_i <= args.max_edit_dis; w_i++){
        correction_process(low_count_reads, w_i);
    }
    
}