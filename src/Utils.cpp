/**
 * @Author: Pengyao Ping
 * @Date:   2023-02-07 12:51:41
 * @Last Modified by:   Pengyao Ping
 * @Last Modified time: 2023-03-05 12:47:18
 */

#include "Utils.hpp"
#include "LoggingLevels.hpp"
#include <iostream>
#include <fstream>
#include <ctime>
#include <iostream>
#include <unordered_map>
#include <set>
#include <map>
#include <thread>
#include <unistd.h>
#include <functional>
#include <chrono>
#include <iomanip>
#include <sstream>
// #include <pbbam/DataSet.h>

int num_cores_to_use;  // Define the global variable

#define reads2graph_VERSION "1.0.0"
#define last_update_date "09.05.2024"

using namespace std;

Utils& Utils::getInstance() {
    static Utils instance;
    return instance;
}

void Utils::logger(int log_level, const std::string& message){
    if (log_level < LOG_LEVEL) return;

    // Get the current time and format it as a string
    time_t rawTime;
    time(&rawTime);
    struct tm * timeInfo = localtime(&rawTime);
    char timeString[20];
    strftime(timeString, 20, "%Y-%m-%d %H:%M:%S", timeInfo);

    std::string log_color;
    std::string log_prefix;
    switch (log_level) {
    case LOG_LEVEL_DEBUG:
        log_color = COLOR_DEBUG;
        log_prefix = "DEBUG: ";
        break;
    case LOG_LEVEL_INFO:
        log_color = COLOR_INFO;
        log_prefix = "INFO: ";
        break;
    case LOG_LEVEL_WARNING:
        log_color = COLOR_WARNING;
        log_prefix = "WARNING: ";
        break;
    case LOG_LEVEL_ERROR:
        log_color = COLOR_ERROR;
        log_prefix = "ERROR: ";
        break;
    case LOG_LEVEL_CRITICAL:
        log_color = COLOR_CRITICAL;
        log_prefix = "CRITICAL ERROR: ";
        break;
    default:
        log_color = COLOR_RESET;
        break;
    }
    logFile << timeString << ": " << log_prefix << message << endl;
    std::cout << timeString << ": " << log_color << log_prefix << COLOR_RESET << message << std::endl;
    // if(logFile.is_open()){
    //     logFile << timeString << ": " << log_prefix << message << endl;
    //     logFile.close();
    // }else {
    //     std::cerr << "\033[1;31mERROR: reads2graph log file opened failed.\033[0m" << std::endl;
    // }
}

Utils::Utils() {
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm now_tm = *std::localtime(&now_c);

    std::ostringstream oss;
    oss << std::put_time(&now_tm, "_%Y%m%d_%H%M%S.log");

    // Get the current time and format it as a string
    time_t rawTime;
    time(&rawTime);
    struct tm * timeInfo = localtime(&rawTime);
    char timeString[20];
    strftime(timeString, 20, "%Y-%m-%d %H:%M:%S", timeInfo);

    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        std::string filePath = std::string(cwd) + "/reads2graph" + oss.str();
        logFile.open(filePath, std::ios::app);
    } else {
        std::cerr << timeString << ": "<< "Error: unable to get current working directory." << std::endl;
    }
}

void Utils::umi_parser(sharg::parser & parser, umi_arguments & args)
{

    // parser.add_option(args.repository, sharg::config{.description = "The repository name to pull from."});
    // parser.add_option(args.branch, sharg::config{.description = "The branch name to pull from."});
 
}

void Utils::read_parser(sharg::parser & parser, read_arguments & args)
{
    // parser.add_option(args.repository, sharg::config{.description = "The repository name to pull from."});
    // parser.add_option(args.branch, sharg::config{.description = "The branch name to pull from."});
 
}

// Define a function to control all the parameters via command line
void Utils::graph_parser(sharg::parser & parser, graph_arguments & args)
{
    parser.info.author = "Pengyao Ping";
    parser.info.short_description = "Error correction and deduplication for a short-read set.";
    parser.info.description.push_back("Accurate and fast error correction and deduplication on short-read sequencing data");
    parser.info.version = reads2graph_VERSION;
    parser.info.date = last_update_date;
    parser.info.short_copyright = "BSD 3-Clause License";
    parser.info.long_copyright = "BSD 3-Clause License\n"
        "Copyright (c) 2024, Pengyao Ping\n"
        "All rights reserved.\n"
        "\n"
        "Redistribution and use in source and binary forms, with or without\n"
        "modification, are permitted provided that the following conditions are met:\n"
        "\n"
        "  1. Redistributions of source code must retain the above copyright notice,\n"
        "     this list of conditions and the following disclaimer.\n"
        "  2. Redistributions in binary form must reproduce the above copyright notice,\n"
        "     this list of conditions and the following disclaimer in the documentation\n"
        "     and/or other materials provided with the distribution.\n"
        "  3. Neither the name of Pengyao Ping nor the names of its contributors may be\n"
        "     used to endorse or promote products derived from this software without\n" 
        "     specific prior written permission.\n"
        "\n"
        "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND \n"
        "ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED \n"
        "WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.\n"
        "IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, \n"
        "INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, \n" 
        "BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,\n"
        " DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY\n"
        "OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE\n"
        "OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED\n"
        "OF THE POSSIBILITY OF SUCH DAMAGE.";

    parser.add_option(args.input_data, 
                        sharg::config{.short_id = 'i', 
                                      .long_id = "input_data", 
                                      .description = "Please provide a fasta/fastq/ data file."});

    // parser.add_option(args.chunk_size,
    //                   sharg::config{.long_id = "chunk_size",
    //                                 .description = "Reading chunk_size records at a time."});

    parser.add_option(args.output_dir,
                      sharg::config{.short_id = 'o',
                                    .long_id = "output_dir",
                                    .description = "The directory for outputs."});

    parser.add_option(args.read_length,
                      sharg::config{.short_id = 'r',
                                    .long_id = "read_length",
                                    .description = "No need to input this parameter, reads2graph will calculate the minimum read length."});

    parser.add_option(args.k_size,
                      sharg::config{.short_id = 'k',
                                    .long_id = "k_size",
                                    .description = "The size for minimiser."});

    parser.add_option(args.window_number,
                      sharg::config{.short_id = 'w',
                                    .long_id = "window_number",
                                    .description = "The window number for minimiser."});

    // parser.add_option(args.window_size,
    //                   sharg::config{.short_id = 'w',
    //                                 .long_id = "window_size",
    //                                 .description = "The window size for minimiser."});

    parser.add_option(args.max_edit_dis,
                      sharg::config{.short_id = 'x',
                                    .long_id = "max_edit_dis",
                                    .description = "The maximum edit distance for constructing edges between reads"});

    parser.add_option(args.min_edit_dis,
                      sharg::config{.short_id = 'n',
                                    .long_id = "min_edit_dis",
                                    .description = "The minimum edit distance for constructing edges between reads."});

    parser.add_option(args.num_process,
                      sharg::config{.short_id = 'p',
                                    .long_id = "num_process",
                                    .description = "The number of expected processes."});

    // parser.add_option(args.graph_filename,
    //                   sharg::config{.short_id = 'g',
    //                                 .long_id = "graph_filename",
    //                                 .description = "The file name of the constructed graph."});
    // Brute Force
    parser.add_option(args.pair_wise,
                      sharg::config{
                                    .long_id = "pair_wise",
                                    .description = "Brute Force calcualte the pairwise edit distance."});

    // parser.add_option(args.bin_size_min,
    //                   sharg::config{.long_id = "bin_size_min",
    //                                 .description = "The smaller threshold used to group buckets of different sizes."});

    parser.add_option(args.bin_size_max,
                      sharg::config{.long_id = "bin_size_max",
                                    .description = "The larger threshold used to group buckets of different sizes."});

    parser.add_option(args.omh_k,
                      sharg::config{.long_id = "omh_k",
                                    .description = "K-mer size used in order min hashing."});

    parser.add_option(args.omh_times,
                      sharg::config{.long_id = "omh_times",
                                    .description = "The number of times to perform permutation in order min hashing."});

    parser.add_option(args.omh_seed,
                      sharg::config{.long_id = "omh_seed",
                                    .description = "The seed to generate a series of seeds for OMH bucketing."});

    parser.add_option(args.omh_flag,
                      sharg::config{.long_id = "omh_flag",
                                    .description = "Do not set this flag by yourself. When the permutation_times larger than the number of k-mer candidates and the kmer size are the same one, bucketing the reads using each kmer candidate."});

    parser.add_option(args.bad_kmer_ratio,
                      sharg::config{.long_id = "bad_kmer_ratio",
                                    .description = "The maximum ratio of bad k-mers out of total number of kmers in a window of a read."});

    parser.add_option(args.probability,
                      sharg::config{.long_id = "probability",
                                    .description = "The expected probability P for grouping two similar reads into same bucket by at least one minimiser that does not include the different bases"});

    parser.add_option(args.visit_depth,
                      sharg::config{.long_id = "visit_depth",
                                    .description = "The maximum distance of nodes from the give node for updating more potential edges."});

    parser.add_option(args.save_graph,
                      sharg::config{.long_id = "save_graph",
                                    .description = "If ture, reads2graph will save graph to file in graphviz dot format."});
    // parser.add_option(args.minimizer_omh,
    //                   sharg::config{.long_id = "minimizer_omh",
    //                                 .description = "If ture, reads2graph employs minimizer bucketing first and then OMH bucketing; otherwise, OMH first then minimizer."});

    // parser.add_option(args.omh_k_step_size,
    //                   sharg::config{.long_id = "omh_k_step_size",
    //                                 .description = "The step size for varied k from the estimated better k for OMH bucketing."});
    // parser.add_option(args.sampling_rate,
    //                   sharg::config{.long_id = "sampling_rate",
    //                                 .description = "Sampling rate for estimating the kmer size."}); 
}


