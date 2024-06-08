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
// #include <gzstream.h>

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
    // void compress_file(const std::string& inputFileName, const std::string& outputFileName);

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

    // Get the stem (file name without extension) and extension of the input file
    std::string stem = input_path.stem().string();
    std::string extension = input_path.extension().string();

    // Check if the input file extension ends with ".gz"
    bool is_gzipped = extension.ends_with(".gz");

    // If the input file is gzipped, remove the ".gz" extension
    if (is_gzipped) {
        extension = input_path.stem().extension().string();
    }

    // Construct the output filename
    std::string output_filename = stem + output_type + extension;

    // Join the output directory path with the output filename
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
    std::string err_correction_flag = ".error_corrected";
    auto err_correction_path_file = get_output_filename(err_correction_flag);
    std::ofstream error_correction_stream(err_correction_path_file);

    std::string deduplication_flag = ".deduplicated";
    auto deduplication_path_file = get_output_filename(deduplication_flag);
    std::ofstream deduplication_stream(get_output_filename(deduplication_flag));

    auto is_fastq = is_fastq_by_extension(args.input_data);

    // Define the type of sequence file output
    if (is_fastq) {
        seqan3::sequence_file_output error_correction_output{error_correction_stream, seqan3::format_fastq{}};

        seqan3::sequence_file_output deduplication_output{deduplication_stream, seqan3::format_fasta{}};

        for (const auto &vertex : boost::make_iterator_range(vertices(graph))) {
            const auto &props = graph[vertex];
            if (args.error_correction){
                for (const auto &id : props.ids) {
                    error_correction_output.emplace_back(props.read, id + " " + id2description[id], id2quality[id]);
                }                
            }
            if (args.deduplication){
                std::string concatenated_ids = std::accumulate(props.ids.begin(), props.ids.end(), std::string{}, [](const std::string& acc, const std::string& id) {
                    return acc.empty() ? id : acc + " " + id;
                });
                deduplication_output.emplace_back(props.read, props.ids[0] + " " + id2description[props.ids[0]] + " " + concatenated_ids);
            }
        }

    } else {
        seqan3::sequence_file_output error_correction_output{error_correction_stream, seqan3::format_fasta{}};
        seqan3::sequence_file_output deduplication_output{deduplication_stream, seqan3::format_fasta{}};
        for (const auto &vertex : boost::make_iterator_range(vertices(graph))) {
            const auto &props = graph[vertex];
            if (args.error_correction){
                for (const auto &id : props.ids) {
                    error_correction_output.emplace_back(props.read, id + " " + id2description[id]);
                }
            }
            if (args.deduplication){
                std::string concatenated_ids = std::accumulate(props.ids.begin(), props.ids.end(), std::string{}, [](const std::string& acc, const std::string& id) {
                    return acc.empty() ? id : acc + " " + id;
                });
                deduplication_output.emplace_back(props.read, props.ids[0] + " " + id2description[props.ids[0]] + " " + concatenated_ids);
            }            
        } 
    }
    // if (args.compression) {
    //     if (args.error_correction){
    //         compress_file(err_correction_path_file, err_correction_path_file + ".gz");
    //     }
    //     if (args.deduplication) {
    //         compress_file(deduplication_path_file, deduplication_path_file + ".gz");
    //     }
    // }

}

// template<typename ArgsType>
// void OutputWriter<ArgsType>::compress_file(const std::string& inputFileName, const std::string& outputFileName) {
//     // Open input file
//     std::ifstream inputFile(inputFileName);
//     if (!inputFile.is_open()) {
//         std::cerr << "Error: Unable to open input file: " << inputFileName << std::endl;
//         return;
//     }

//     // Open output file
//     ogzstream outputFile(outputFileName.c_str());
//     if (!outputFile.is_open()) {
//         std::cerr << "Error: Unable to open output file: " << outputFileName << std::endl;
//         return;
//     }

//     // Read input file line by line and write to output gzipped file
//     std::string line;
//     while (std::getline(inputFile, line)) {
//         outputFile << line << '\n';
//     }

//     // Close files
//     inputFile.close();
//     outputFile.close();
//     std::filesystem::remove(inputFile);
// }

#endif