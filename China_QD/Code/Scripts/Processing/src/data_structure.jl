
struct ParsedData
    channel::Vector{UInt32}
    delays::Vector{Float64}
    arrivals::Vector{Float64}
end

struct BinnedData
    intensity::Vector{UInt32}
    index_min::Vector{UInt32}
    index_max::Vector{UInt32}
end

struct FLIDData
    tau_F::Vector{Float64}
    tau_D::Vector{Float64}
    A2::Vector{Float64}
    bkg::Vector{Float64}
    I::Vector{Float64}
end