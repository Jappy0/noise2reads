// MinimizerGenerator.tpp

template<typename ArgsType>
MinimizerGenerator<ArgsType>::MinimizerGenerator(ArgsType args) : args(args) {}

template<typename ArgsType>
std::unordered_map<std::uint64_t, std::vector<std::vector<seqan3::dna5>>> MinimizerGenerator<ArgsType>::minimizer2reads_main(std::vector<std::vector<seqan3::dna5>> unique_reads,std::tuple<unsigned, unsigned, unsigned, double> betterParams)
{   
    // auto better_n = std::get<0>(betterParams);
    auto better_ww = std::get<1>(betterParams);
    auto better_kk = std::get<2>(betterParams);
    // auto prob = std::get<3>(betterParams);
    // auto best_w = round(static_cast<double>(args.read_length) / best_n);   
    auto better_k = static_cast<uint8_t>(better_kk);
    auto better_w = static_cast<uint8_t>(better_ww);

    // int available_cores = omp_get_max_threads();
    // auto num_cores_to_use = std::min(std::max(args.num_process, 1), available_cores);
    // omp_set_num_threads(num_cores_to_use);

    // #pragma omp parallel for
    // for (size_t i = 0; i < read2count.size(); ++i) {
    //     auto it = std::next(read2count.begin(), i);
    //     const auto& [read, count] = *it;
    // #pragma omp parallel for
    #pragma omp parallel for num_threads(args.num_process) schedule(static)
    for (auto const & read : unique_reads){
        // auto minimisers = read | seqan3::views::kmer_hash(seqan3::ungapped{args.k_size}) | seqan3::views::minimiser(args.window_size - args.k_size + 1);
        auto minimisers = read | seqan3::views::kmer_hash(seqan3::ungapped{better_k}) | seqan3::views::minimiser(better_w - better_k + 1);   

        // // Iterate over minimisers and group reads
        for (auto const &minimiser : minimisers) {
            std::uint64_t converted_minimiser = static_cast<std::uint64_t>(minimiser);
            #pragma omp critical
            {
                minimiser_to_reads[converted_minimiser].push_back(read);
            }
        }
    }

    Utils::getInstance().logger(LOG_LEVEL_DEBUG, boost::str(boost::format("Size of minimiser_to_reads: %1%!") % minimiser_to_reads.size())); 
    return minimiser_to_reads;     
}

template<typename ArgsType>
int MinimizerGenerator<ArgsType>::kSize(int L, double p) {
    return ceil((p*(1+L))/(1+p));
}

template<typename ArgsType>
double MinimizerGenerator<ArgsType>::proba(unsigned L, unsigned k) {
    double p;
    p = (static_cast<double>(k))/(L-k+1);
    return p;
}

template<typename ArgsType>
std::tuple<unsigned, unsigned, unsigned, double> MinimizerGenerator<ArgsType>::possibleBetterParameters() {
    unsigned betterK;
    unsigned betterN;
    unsigned betterW;
    double p=0;
    // if (args.read_length >= 8 && args.read_length < 16){
    //     betterK = 3;
    //     betterN = 2;
    //     // betterW = args.read_length;
    //     betterW = round(args.read_length/betterN);
    // } else if (args.read_length >= 16 && args.read_length < 50){
    //     betterK = 4;
    //     betterN = 2;
    //     betterW = round(args.read_length/betterN);
    // } else if (args.read_length >= 50 && args.read_length <= 300) {
    //     // if (args.max_edit_dis == 1 || args.max_edit_dis == 2){
    //     //     betterN = 3;
    //     // } else {
    //     //     betterN = ceil((static_cast<double>(args.max_edit_dis))/2)+1;
    //     // }
    //     betterN = 3;
    //     betterW = round(args.read_length/betterN);
    //     betterK = kSize(betterW, args.bad_kmer_ratio);
    // } 
    if (args.read_length >= 6 && args.read_length < 10){
        betterK = args.k_size;
        betterN = args.window_number;
        betterW = args.read_length;
    } else if (args.read_length >= 10 && args.read_length < 16){
        betterK = args.k_size;
        betterN = args.window_number;
        // betterW = args.read_length;
        betterW = round(args.read_length/betterN);
    } else if (args.read_length >= 16 && args.read_length < 50){
        betterK = args.k_size;
        betterN = args.window_number;
        betterW = round(args.read_length/betterN);
    } else if (args.read_length >= 50 && args.read_length <= 300) {
        // if (args.max_edit_dis == 1 || args.max_edit_dis == 2){
        //     betterN = 3;
        // } else {
        //     betterN = ceil((static_cast<double>(args.max_edit_dis))/2)+1;
        // }
        betterN = args.window_number;
        betterW = round(args.read_length/betterN);
        betterK = kSize(betterW, args.bad_kmer_ratio);
        if (betterK < 4){
            Utils::getInstance().logger(LOG_LEVEL_WARNING, boost::str(boost::format("Estimated k=%1% has been changed to 4.") % betterK));
            betterK = 4;
        } else 
        if (betterK >= 28) {
            Utils::getInstance().logger(LOG_LEVEL_WARNING, boost::str(boost::format("Estimated k=%1% has been changed to 27 as the maximum size of unggaped shape is stricted by 28 in Seqan3.") % betterK));  
            betterK = 27;             
        }
    } 
    // if (betterK < 4){
    //     Utils::getInstance().logger(LOG_LEVEL_WARNING, std::format("Better k {} has been changed to 4.", betterK));
    //     betterK = 4;
    // } else 
    // if (betterK >= 28) {
    //     Utils::getInstance().logger(LOG_LEVEL_WARNING, std::format("Better k {} has been changed to 27 as the maximum size of unggaped shape is stricted by 28 in Seqan3.", betterK));  
    //     betterK = 27;             
    // }
    if (args.read_length >= 50 && args.read_length <= 300){
        p = 1 - std::pow(proba(betterW, betterK), betterN);
        Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("Estimated number of windows: %1%, Estimated window size: %2%, Estimated K: %3% and the probability: %4%.") % betterN % betterW % betterK % p)); 
    } else {
        Utils::getInstance().logger(LOG_LEVEL_INFO, boost::str(boost::format("Number of windows: %1%, Window size: %2%, K size: %3%.") % betterN % betterW % betterK));        
    }
    auto number_kmer = betterW - betterK + 1;
    if ((number_kmer) < 3 ){
        Utils::getInstance().logger(LOG_LEVEL_WARNING, boost::str(boost::format("only %1% kmers setted for minimizer selection.") % number_kmer));
    }

    return std::make_tuple(betterN, betterW, betterK, p);   
}