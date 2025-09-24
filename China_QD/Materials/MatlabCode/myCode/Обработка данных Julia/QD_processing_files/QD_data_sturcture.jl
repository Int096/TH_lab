struct ParsedData
    channel::Vector
    delays::Vector
    arrivals::Vector
end

struct BinnedData
    intensity::Vector{UInt32}
    index_min::Vector{UInt32}
    intex_max::Vector{UInt32}
end

struct CorrData
    g2::Vector
    times::Vector 
    CoinCounts::Vector 
    new_t::Vector
    new_p::Vector
    new_N_kor
    experiment_time 
    new_delta_t 
    n_m  
end