/*
 * @Author: Pengyao PING
 * @Date: 2021-07-19 17:24:08
 * @LastEditors: Pengyao PING
 * @LastEditTime: 2021-08-11 11:36:10
 * @Email: Pengyao.Ping@student.uts.edu.au
 * @Description: 
 */

#include "Utils.hpp"
#include "Parser.hpp"
#include "LoggingLevels.hpp"
#include "GraphConstructor.hpp"
#include "Reader.hpp"
#include "MinimizerGenerator.hpp"
#include "OMH.hpp"
#include "ReadCorrection.hpp"
#include "OutputWriter.hpp"

#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <getopt.h>
#include <fstream>
#include <filesystem>
#include <unistd.h>
#include <seqan3/core/debug_stream.hpp> // for debug_stream

#include <array>  // std::array
#include <string> // std::string
#include <vector> // std::vector
 
#include <seqan3/alphabet/all.hpp>
#include <seqan3/alphabet/views/all.hpp> // optional: use views to convert the input string to a dna5 sequence
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/utility/views/all.hpp> // optional: use views to convert the input string to a dna5 sequence
#include <seqan3/io/sequence_file/all.hpp>
// #include <seqan3/std/filesystem>
#include <boost/format.hpp>
#include <omp.h>

using namespace std;
using namespace seqan3::literals;

int main(int argc, char** argv) {
    Utils::getInstance().logger(LOG_LEVEL_INFO,  "Welcome to use noise2reads!");
    /////////
    std::string concatenatedArgs;
    for (int i = 0; i < argc; ++i) {
        concatenatedArgs += argv[i];
        if (i < argc - 1) {
            concatenatedArgs += " ";
        }
    }
    Utils::getInstance().logger(LOG_LEVEL_INFO, concatenatedArgs);
    //////////////////
    // Declare and define a global variable for available cores
    int available_cores = omp_get_max_threads();
    Utils::getInstance().logger(LOG_LEVEL_DEBUG, boost::str(boost::format("The maximum number of threads available: %1% ") % available_cores));
    /////////////////////////////////////////////////////////////////
    sharg::parser top_level_parser{"noise2reads", argc, argv, sharg::update_notifications::off, {"umi", "read", "graph"}};

    // Add information and flags, but no (positional) options to your top-level parser.
    // Because of ambiguity, we do not allow any (positional) options for the top-level parser.
    // top_level_parser.info.description.push_back("Error correction and deduplication on short-read sequencing data.");
    // bool flag{false};
    // top_level_parser.add_flag(flag, sharg::config{.short_id = 'f', .long_id = "flag", .description = "some flag"});

    try
    {
        top_level_parser.parse(); // trigger command line parsing
    }
    catch (sharg::parser_error const & ext) // catch user errors
    {
        std::cerr << "[Error] " << ext.what() << "\n"; // customise your error message
        return -1;
    }

    sharg::parser & sub_parser = top_level_parser.get_sub_parser();

    /////////////////////////////////////////////////////////////////////////////////////////////

    if (sub_parser.info.app_name == std::string_view{"noise2reads-umi"} || sub_parser.info.app_name == std::string_view{"noise2reads-read"}){
        using Args_Type = std::conditional_t<std::is_same_v<decltype(sub_parser.info.app_name), std::string_view>, umi_arguments, read_arguments>;
        Args_Type args;
        Parser<Args_Type>().read_parser(sub_parser, args);
        try
        {
            sub_parser.parse(); // trigger command line parsing
        }
        catch (sharg::parser_error const & ext) // catch user errors
        {
            std::cerr << "[Invalid Options] " << ext.what() << "\n"; // customise your error message
            return -1;
        }

        args.num_process = std::min(std::max(args.num_process, 1), available_cores);
        Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of threads: %1% ") % args.num_process));
        //////////////////////////////////////////////////////////////////
        // construct an undirected edit-distance-based read graph
        Reader<Args_Type> reader(args);
        auto [unique_reads, read2count, read2ids, min_read_length] = reader.get_unique_reads_counts();
        args.read_length = min_read_length;
        auto total_uniq_num = unique_reads.size();
        Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of unique reads: %1%, minimum read length: %2%.") % total_uniq_num % min_read_length));
        //////////////////////////////////////////////////////////////
        Graph nt_ed_graph;
        GraphConstructor<Args_Type> graph_constructor(read2count, read2ids, args);
        if (args.pair_wise) {
            nt_ed_graph = graph_constructor.construt_graph_via_pairwise_comparison(unique_reads);
            if (args.save_graph){
                graph_constructor.save_graph();
            }
        } else {
            // minimizer grouping first and then omh
            MinimizerGenerator<Args_Type> minimizer_generator(args);
            auto betterParams = minimizer_generator.possibleBetterParameters();
            auto hash2reads = minimizer_generator.minimizer2reads_main(unique_reads, betterParams);  

            nt_ed_graph = graph_constructor.construct_graph(hash2reads);
            if (args.save_graph){
                graph_constructor.save_graph();
            }
        }

        ReadCorrection read_correction(nt_ed_graph, args, graph_constructor.get_read2vertex(), graph_constructor.get_vertex2read());
        read_correction.correction_main(read2count);
        Graph update_ed_graph = read_correction.get_graph();

        auto id2quality = reader.get_id2quality();
        auto id2description = reader.get_id2description();
        OutputWriter<Args_Type>(args).write_graph2dataset(std::move(update_ed_graph), id2quality, id2description);
        //////////////////////////////////////////////////////////////            
    } else if (sub_parser.info.app_name == std::string_view{"noise2reads-graph"}){
        using Args_Type = graph_arguments;
        Args_Type args;
        Parser<Args_Type>().graph_parser(sub_parser, args);
        try
        {
            sub_parser.parse(); // trigger command line parsing
        }
        catch (sharg::parser_error const & ext) // catch user errors
        {
            std::cerr << "[Invalid Options] " << ext.what() << "\n"; // customise your error message
            return -1;
        }
        ////////////////////////////////////////////////////////////////////////////
        // Ensure the user-specified number of cores is within a valid range
        args.num_process = std::min(std::max(args.num_process, 1), available_cores);
        Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of threads: %1% ") % args.num_process));

        ////////////////////////////////////////////////////////////
        auto [unique_reads, read2count, read2ids, min_read_length] = Reader<Args_Type>(args).get_unique_reads_counts();
        args.read_length = min_read_length;
        auto total_uniq_num = unique_reads.size();
        Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("The number of unique reads: %1%, minimum read length: %2%.") % total_uniq_num % min_read_length));
        //////////////////////////////////////////////////////////////
        GraphConstructor<Args_Type> graph_constructor(read2count, read2ids, args);
        if (args.pair_wise) {
            graph_constructor.construt_graph_via_pairwise_comparison(unique_reads);
            if (args.save_graph){
                graph_constructor.save_graph();
            }
        } else {
            // minimizer grouping first and then omh
            MinimizerGenerator<Args_Type> minimizer_generator(args);
            auto betterParams = minimizer_generator.possibleBetterParameters();

            auto hash2reads = minimizer_generator.minimizer2reads_main(unique_reads, betterParams);  

            graph_constructor.construct_graph(hash2reads);
            if (args.save_graph){
                graph_constructor.save_graph();
            }
        }
        ////////////////////////////////////////////////////////////////////////////
    } else {
        std::cout << "Unhandled subparser named " << sub_parser.info.app_name << '\n';
    }
    return 0;
}
