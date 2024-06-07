#ifndef READWRITE_HPP
#define READWRITE_HPP

#include "Utils.hpp"
#include "LoggingLevels.hpp"
#include <array>  // std::array
#include <string> // std::string
#include <vector> // std::vector
 
#include <seqan3/alphabet/all.hpp>
#include <seqan3/alphabet/views/all.hpp> // optional: use views to convert the input string to a dna5 sequence
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/utility/views/all.hpp> // optional: use views to convert the input string to a dna5 sequence
#include <seqan3/io/sequence_file/all.hpp>

using namespace std;

template<typename ArgsType>
class ReadWrite{
    public:
        ReadWrite(ArgsType args);//constructor
        ~ReadWrite(); //deconstructor
        std::tuple<std::vector<std::vector<seqan3::dna5>>, std::map<std::vector<seqan3::dna5>, uint32_t>, std::map<std::vector<seqan3::dna5>, std::vector<std::string>>, unsigned> get_unique_reads_counts();
    private:
        // graph_arguments args;
        ArgsType args;
};

#include "../src/ReadWrite.tpp"

#endif