// ReadCorrection.h
#ifndef ReadCorrection_HPP
#define ReadCorrection_HPP

#include "GraphConstructor.hpp"
#include "Parser.hpp"

#include <boost/property_map/property_map.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <boost/format.hpp>

using namespace std;
using namespace boost;
using namespace seqan3::literals;

class ReadCorrection
{
public:
    ReadCorrection(Graph graph, read_arguments args, std::unordered_map<std::vector<seqan3::dna5>, Vertex, std::hash<std::vector<seqan3::dna5>>> read2vertex, std::unordered_map<Vertex, std::vector<seqan3::dna5>, std::hash<Vertex>> vertex2read);
    std::vector<std::vector<seqan3::dna5>> get_low_frequency_reads(const std::map<std::vector<seqan3::dna5>, uint32_t>& read2count);
    bool is_isolated(Vertex v);
    std::vector<Vertex> get_connected_nodes_with_weight(Vertex v, int w);
    void correction_process(const std::vector<std::vector<seqan3::dna5>>& low_count_reads, int w);
    void correction_main(const std::map<std::vector<seqan3::dna5>, uint32_t>& read2count);
    void correction_isolates(const std::map<std::vector<seqan3::dna5>, uint32_t>& read2count);

    void update_graph(std::vector<std::vector<seqan3::dna5>> high_freq_reads, int w);
    void visitNeighborsWithThreshold(const Graph& g, Vertex node, int distance_threshold, int current_distance, std::vector<Vertex>& indirect_neighbors, std::vector<bool>& visited);
    std::vector<Vertex> visitNeighborsOfNeighborsWithThreshold(const Graph& g, Vertex node, int distance_threshold);
    std::pair<std::vector<std::vector<seqan3::dna5>>, std::vector<std::vector<seqan3::dna5>>> get_high_low_freq_reads();
    void insert_edge(std::vector<seqan3::dna5> read1, std::vector<seqan3::dna5> read2, int edit_dis);

    Graph get_graph() const {
        return graph_;
    }

private:
    // std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> key2reads_;
    // std::map<std::vector<seqan3::dna5>, uint32_t> read2count_;
    Graph graph_;
    read_arguments args;
    std::unordered_map<std::vector<seqan3::dna5>, Vertex, std::hash<std::vector<seqan3::dna5>>> read2vertex_;
    std::unordered_map<Vertex, std::vector<seqan3::dna5>, std::hash<Vertex>> vertex2read_;
};

#endif // GraphConstructor_H