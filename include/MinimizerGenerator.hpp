// MinimizerGenerator.hpp
#ifndef __MINIMIZERGENERATOR_HPP__
#define __MINIMIZERGENERATOR_HPP__

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/all.hpp>
#include <vector>
#include <set>
#include <seqan3/alphabet/all.hpp>
#include <algorithm>
#include <iostream>
#include <type_traits> // for std::decay_t
#include "Utils.hpp"
#include "LoggingLevels.hpp"

using namespace std;
using namespace seqan3::literals;

class MinimizerGenerator
{
public:
    MinimizerGenerator(graph_arguments args);
    // MinimizerGenerator(std::map<std::vector<seqan3::dna5>, uint32_t> read2count, graph_arguments args);
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> minimizer2reads_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::tuple<unsigned, unsigned, unsigned, double> betterParams);
    double proba(unsigned L, unsigned k);
    int kSize(int L, double p);
    std::tuple<unsigned, unsigned, unsigned, double> possibleBetterParameters();
private:
    // std::vector<std::vector<seqan3::dna5>> unique_reads;
    // std::map<std::vector<seqan3::dna5>, uint32_t> read2count;
    std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> minimiser_to_reads;
    graph_arguments args;
};

#endif