#ifndef OUTPUTWRITER_HPP
#define OUTPUTWRITER_HPP

#include "Reader.hpp"
#include "GraphConstructor.hpp"
#include "ReadCorrection.hpp"

#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <seqan3/alphabet/nucleotide/dna5.hpp>  // Include SeqAn header as needed
#include <seqan3/alphabet/all.hpp>
#include <seqan3/alphabet/views/all.hpp> 
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/utility/views/all.hpp> 
#include <seqan3/io/sequence_file/all.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <boost/format.hpp>
#include <boost/graph/adjacency_list.hpp> 

using namespace std;
using namespace boost;
using namespace seqan3::literals;

template<typename ArgsType>
class OutputWriter {
public:
    OutputWriter(ArgsType args);
    ~OutputWriter(void);
    void write_graph2dataset(Graph graph, std::map<std::string, std::vector<seqan3::phred42>> id2quality, std::map<std::string, std::string> id2description);
    std::string get_output_filename(std::string output_type);
    bool is_fastq_by_extension(const std::string &input_file);

private:
    ArgsType args;
};

template<typename ArgsType>
OutputWriter<ArgsType>::OutputWriter(ArgsType args) : args(args){}

template<typename ArgsType>
OutputWriter<ArgsType>::~OutputWriter(void){}

template<typename ArgsType>
std::string OutputWriter<ArgsType>::get_output_filename(std::string output_type) {
    std::filesystem::path input_path(args.input_data);
    std::filesystem::path output_path(args.output_dir);
    std::string output_filename = input_path.stem().string() + output_type + input_path.extension().string();
    return (output_path / output_filename).string();
}

template<typename ArgsType>
bool OutputWriter<ArgsType>::is_fastq_by_extension(const std::string &input_file) {
    std::string ext = std::filesystem::path(input_file).extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower); // Convert to lower case for case-insensitive comparison

    if (ext == ".gz") {
        std::string stem_ext = std::filesystem::path(input_file).stem().extension().string();
        std::transform(stem_ext.begin(), stem_ext.end(), stem_ext.begin(), ::tolower); // Convert to lower case for case-insensitive comparison
        
        if (stem_ext == ".fastq" || stem_ext == ".fq") {
            return true;
        } else if (stem_ext == ".fasta" || stem_ext == ".fa") {
            return false;
        } else {
            throw std::runtime_error("Unknown file extension before .gz: " + stem_ext);
        }
    } else {
        if (ext == ".fastq" || ext == ".fq") {
            return true;
        } else if (ext == ".fasta" || ext == ".fa") {
            return false;
        } else {
            throw std::runtime_error("Unknown file extension: " + ext);
        }
    }
}

template<typename ArgsType>
void OutputWriter<ArgsType>::write_graph2dataset(Graph graph, std::map<std::string, std::vector<seqan3::phred42>> id2quality, std::map<std::string, std::string> id2description){
    // Define the type of sequence file input
    seqan3::sequence_file_input input{args.input_data};
    std::string out_str = ".correction";
    auto output_correction_file = get_output_filename(out_str);
    std::ofstream output_stream(output_correction_file);

    auto is_fastq = is_fastq_by_extension(args.input_data);

    // Define the type of sequence file output
    if (is_fastq) {
        seqan3::sequence_file_output output{output_stream, seqan3::format_fastq{}};
        for (const auto &vertex : boost::make_iterator_range(vertices(graph))) {
            const auto &props = graph[vertex];
            for (const auto &id : props.ids) {
                output.emplace_back(props.read, id + " " + id2description[id], id2quality[id]);
            }
        }
    } else {
        seqan3::sequence_file_output output{output_stream, seqan3::format_fasta{}};
        for (const auto &vertex : boost::make_iterator_range(vertices(graph))) {
            const auto &props = graph[vertex];
            for (const auto &id : props.ids) {
                output.emplace_back(props.read, id + " " + id2description[id]);
            }
        } 
    }
}

#endif