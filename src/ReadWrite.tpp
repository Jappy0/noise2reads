/**
 * @brief Construct a new Read Write:: Read Write object
 * 
 */
template<typename ArgsType>
ReadWrite<ArgsType>::ReadWrite(ArgsType args) : args(args){}

template<typename ArgsType>
ReadWrite<ArgsType>::~ReadWrite(void){}

template<typename ArgsType>
std::tuple<std::vector<std::vector<seqan3::dna5>>, std::map<std::vector<seqan3::dna5>, uint32_t>, std::map<std::vector<seqan3::dna5>, std::vector<std::string>>, unsigned> ReadWrite<ArgsType>::get_unique_reads_counts(){

    seqan3::sequence_file_input fin{args.input_data};
    // using record_type = decltype(fin)::record_type;
    // std::vector<record_type> records{};
    // std::vector<decltype(fin)::record_type> records;
    // Define a set to store unique reads.
    std::vector<std::vector<seqan3::dna5>> unique_reads;
    // Define a map to store unique reads and their counts
    std::map<std::vector<seqan3::dna5>, uint32_t> read2count;
    std::map<std::vector<seqan3::dna5>, std::vector<std::string>> read2ids;
    unsigned min_read_length =  std::numeric_limits<unsigned>::max();

    for (auto &record : fin)
    {
        // records.push_back(std::move(record));
        std::vector<seqan3::dna5> cur_seq = record.sequence();
        read2ids[cur_seq].push_back(record.id());

        if (!read2count[cur_seq]++) {
            unique_reads.push_back(cur_seq);
            auto cur_read_len = cur_seq.size();
            if (cur_read_len < 6 || cur_read_len > 300){
                Utils::getInstance().logger(LOG_LEVEL_ERROR,  "Read length should not less than 6 and larger than 300.");
            }
            min_read_length = std::min(min_read_length, static_cast<unsigned>(cur_read_len));
        }
    }

    Utils::getInstance().logger(LOG_LEVEL_INFO,  "Loading data done!");
    return {unique_reads, read2count, read2ids, min_read_length};   
}