// OMH.tpp
template<typename ArgsType>
OMH<ArgsType>::OMH(ArgsType args) : args(args) {}

template<typename ArgsType>
unsigned OMH<ArgsType>::omh_k(unsigned L, double p, uint8_t d) {
    // unsigned k = ceil((p*(1+L))/(d+p));
    // unsigned k = ceil(((1-p)*(1+L))/(d+1-p));
    unsigned k;
    if (args.read_length >= 6 && args.read_length < 50){
        auto omh_kmer_n = args.read_length- 2 * args.omh_k + 1;
        if (omh_kmer_n < 3 ){
            Utils::getInstance().logger(LOG_LEVEL_WARNING,  std::format("only {} kmers setted for gOMH selection.", omh_kmer_n));
        }        
        k = args.omh_k;
    } else if (args.read_length >= 50 && args.read_length <= 300){
        k = ceil(((1-p)*(2+L))/(d+2-2*p));
        if (k < 4){
            Utils::getInstance().logger(LOG_LEVEL_WARNING, std::format("Estimated k={} has been changed to 4.", k));
            k = 4;
        } else if (k > 27) {
            Utils::getInstance().logger(LOG_LEVEL_WARNING, std::format("Estimated k={} has been changed to 27.", k)); 
            k = 27;               
        }   
    }
    return k;
}

template<typename ArgsType>
std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> OMH<ArgsType>::omh2reads_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::uint64_t seed, unsigned k){
    #pragma omp parallel for num_threads(args.num_process) schedule(static)
    for (auto const & read : unique_reads){
        auto omh_value = omh_pos(read, k, seed); 
        #pragma omp critical
        {
            omh2reads[omh_value].push_back(read);
        }        
    }
    Utils::getInstance().logger(LOG_LEVEL_DEBUG, "All the unique reads for OMH done.");  
    return omh2reads;            
}

template<typename ArgsType>
std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> OMH<ArgsType>::omh2read_main(std::vector<std::vector<seqan3::dna5>> unique_reads, std::vector<std::pair<std::uint64_t, unsigned>> seeds_k){
    if (args.omh_flag){
        // When the permutation_times larger than the number of k-mer candidates and the kmer size are the same one, bucketing the reads using each kmer candidate
        auto first_pair = seeds_k[0];
        std::uint64_t seed = first_pair.first;
        unsigned k = first_pair.second;
        #pragma omp parallel for num_threads(args.num_process) schedule(static)
        for (auto const & read : unique_reads){  
            auto omh_values = omh_pos2(read, k, seed);
            for (auto const & omh_val : omh_values){
                #pragma omp critical
                {
                    omh2reads[omh_val].push_back(read);
                }                    
            }
        }  
    } else {
        #pragma omp parallel for num_threads(args.num_process) schedule(static)
        for (auto const & read : unique_reads){
            #pragma omp parallel for num_threads(args.num_process) schedule(static)
            for(auto &pair : seeds_k){
                std::uint64_t seed = pair.first;
                unsigned k = pair.second;
                auto omh_value = omh_pos(read, k, seed); 
                #pragma omp critical
                {
                    omh2reads[omh_value].push_back(read);
                }        
            } 
        }        
    }
    Utils::getInstance().logger(LOG_LEVEL_DEBUG, "All the unique reads for OMH done.");  
    return omh2reads;            
}

template<typename ArgsType>
std::string OMH<ArgsType>::getGappedSubstring(const std::string& str, size_t startPos, size_t length) {
    std::string gappedSubstring;

    for (size_t i = startPos, count = 0; i < str.size() && count < length; i=i+2) {
        gappedSubstring += str[i];
        ++count;
    }

    return gappedSubstring;
}

template<typename ArgsType>
uint64_t OMH<ArgsType>::omh_pos(const std::vector<seqan3::dna5>& read, unsigned k, std::uint64_t seed) {
    if(read.size() < 2*k - 1) return {};
    std::vector<std::uint64_t> hash_vec;
    std::unordered_map<std::string, unsigned> occurrences;
    std::uint64_t cur_seed = seed;

    auto seql = read | seqan3::views::to_char;
    string read_str(seql.begin(), seql.end());

    size_t ll = read_str.size();

    // for(size_t i = 0; i < read_str.size() - k + 1; ++i) {
    //////////////////
    // for(size_t i = 0; i < ll; i += k) {
    //     string kmer;
    //     size_t remaining_length = ll - i;

    //     if (remaining_length >= k) {
    //         kmer = read_str.substr(i, k);
    //     } else if (remaining_length >= k/2){
    //         kmer = read_str.substr(i);
    //     } 
    //     occurrences[kmer]++;
    //     boost::hash_combine(cur_seed, kmer);
    //     boost::hash_combine(cur_seed, occurrences[kmer]);
    //     hash_vec.emplace_back(cur_seed);
    //     cur_seed = seed;
    // }

    for(size_t i = 0; i <= ll - 2*k + 1; ++i) {
        string kmer = getGappedSubstring(read_str, i, k);
        occurrences[kmer]++;
        boost::hash_combine(cur_seed, kmer);
        boost::hash_combine(cur_seed, occurrences[kmer]);
        hash_vec.emplace_back(cur_seed);
        cur_seed = seed;
    }

    auto min_hash = std::min_element(hash_vec.begin(), hash_vec.end());    
    return *min_hash;
}

// using each of all the kmers for bucketing

template<typename ArgsType>
std::vector<std::uint64_t> OMH<ArgsType>::omh_pos2(const std::vector<seqan3::dna5>& read, unsigned k, std::uint64_t seed) {
    if(read.size() < 2*k - 1) return {};
    std::vector<std::uint64_t> hash_vec;
    std::unordered_map<std::string, unsigned> occurrences;
    std::uint64_t cur_seed = seed;

    auto seql = read | seqan3::views::to_char;
    string read_str(seql.begin(), seql.end());

    size_t ll = read_str.size();

    for(size_t i = 0; i <= ll - 2*k + 1; ++i) {
        string kmer = getGappedSubstring(read_str, i, k);
        occurrences[kmer]++;
        boost::hash_combine(cur_seed, kmer);
        boost::hash_combine(cur_seed, occurrences[kmer]);
        hash_vec.emplace_back(cur_seed);
        cur_seed = seed;
    }
    return hash_vec;
}