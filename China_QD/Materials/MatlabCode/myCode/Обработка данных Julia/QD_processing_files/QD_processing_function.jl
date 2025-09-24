using Statistics 
using StatsBase

include("QD_data_sturcture.jl")
include("QD_helper_function.jl")

function ParsingData(source_data::Vector{UInt32})
    # Разбор 32-битных слов на компоненты
    control_bit = UInt8.((source_data .>> 31) .& 0x01)
    channel     = UInt8.((source_data .>> 25) .& 0x3F)
    delays      = UInt64.((source_data .>> 10) .& 0x7FFF)
    impulse_num = UInt64.(source_data .& 0x3FF)

    # Создание массива для корректировки переполнений
    overflow_idxs = findall(==(1), control_bit)
    overflow_diff = Statistics.diff(overflow_idxs) .- 1
    overflow_vec  = collect(0 : length(overflow_diff)-1)

    overflow = Vector{UInt64}(undef, sum(overflow_diff))
    i = 1
    @inbounds for j in eachindex(overflow_vec)
        value   = overflow_vec[j]
        repeats = overflow_diff[j]
        for k in 1 : repeats
            overflow[i] = value
            i += 1
        end
    end
    overflow .*= 1023

    arrivals = impulse_num

    # Удаление битов с control_bit = 1
    mask = ones(Bool, length(control_bit))
    mask[overflow_idxs] .= false

    channel  = channel[mask]
    delays   = delays[mask]
    arrivals = arrivals[mask]
    impulse_num = impulse_num[mask]

    # Корректировка битов переполнения и задержек
    impulse_num .+= overflow
    delays .*= 64

    # Формирование arrivals массива
    time_between = maximum(delays)
    arrivals = (delays .+ impulse_num .* UInt64(time_between)) .* 1e-12
    arrivals .-= arrivals[1]
    delays = delays .* 1e-12

    # Сортировка
    perm = sortperm(arrivals)
    delays = delays[perm]
    channel = channel[perm]
    
    return ParsedData(channel, delays, arrivals)
end

function BinningData(arrivals::Vector{Float64}, bin_size::Float64)
    nums_of_bin = UInt32(arrivals[end] ÷ bin_size)
    nums_of_arrivals = length(arrivals)

    intensity = Vector{UInt32}(undef, nums_of_bin)
    index_min = Vector{UInt32}(undef, nums_of_bin)
    index_max = Vector{UInt32}(undef, nums_of_bin)

    count = 0
    for j in 1:nums_of_bin
        if j == 1
            index_min[j] = 1
            count = 1
        else
            index_min[j] = index_max[j-1] + 1
            count = index_min[j]
        end

        while arrivals[count] < j*bin_size
            count += 1
            if count + 1 > nums_of_arrivals
                break
            end
        end

        if count == index_min[j]
            intensity[j] = 0
            index_max[j] = index_min[j] - 1
        else
            index_max[j] = count - 1
            intensity[j] = index_max[j] - index_min[j] + 1
        end
    end

    return BinnedData(intensity, index_min, index_max)
end

function CalculateCorrelateFunction(parsed_data::ParsedData, intensity::Vector{UInt32})
    arrivals = parsed_data.arrivals
    delays   = parsed_data.delays
    channel  = parsed_data.channel

    ch1 = arrivals[channel .== 0]
    ch2 = arrivals[channel .== 1]

    nums_of_points::Int32 = 10
    time_between::Float64 = maximum(delays)
    end_coef::Float64     = 0.75

    g2, times, CoinCounts = GetCorrelation(ch1, ch2, nums_of_points, time_between, end_coef)

    t = times'
    p = (g2 .- 1)'

    experiment_time = max(maximum(ch1), maximum(ch2)) - min(minimum(ch1), minimum(ch2))

    h = fit(Histogram, intensity, minimum(intensity):maximum(intensity))
    Y_exp = h.weights / sum(h.weights)
    X_exp = collect(h.edges[1])[1:end-1]
    Y_exp = Y_exp / sum(Y_exp)
    N_med = Y_exp * X_exp'

        new_bins = vcat(
        exp10.(range(-4.5, -4, length=2)),
        exp10.(range(-3.5, -3, length=2)),
        exp10.(range(-2.6667, -2, length=3)),
        exp10.(range(-1.6667, -1, length=3)),
        exp10.(range(-0.8, 0, length=5)),
        exp10.(range(0.1, 1, length=10)),
        exp10.(range(1.0667, 2, length=10)),
        exp10.(range(2.1, log10(maximum(t)), length=10))
    )
    num_new_bins = length(new_bins)

    new_p = zeros(Float64, num_new_bins)
    new_t = zeros(Float64, num_new_bins)

    idx_photons_in_bins = (t .<= new_bins[1]) .& (t .> 0)
    new_p[1] = mean(p[idx_photons_in_bins])
    new_t[1] = mean(t[idx_photons_in_bins])

    for ii in 2:num_new_bins
        idx_photons_in_bins = (t .<= new_bins[ii]) .& (t .> new_bins[ii-1])
        new_p[ii] = mean(p[idx_photons_in_bins])
        new_t[ii] = mean(t[idx_photons_in_bins])
    end

    new_delta_t = diff(vcat([0.0], new_t))
    n_m = length(ch2) / experiment_time
    new_N_kor = (new_p .+ 1) .* new_delta_t * n_m^2 * experiment_time

    return CorrData(g2, times, CoinCounts, new_t, new_p, new_N_kor, experiment_time, new_delta_t, n_m)
end