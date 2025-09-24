include("data_structure.jl")

using Statistics, StatsBase, JLD2

function MakeDirs(folder::String, res_fol::String, file::String)
    path_to_save_folder = joinpath(res_fol, folder, file)
    ispath(path_to_save_folder) || mkdir(path_to_save_folder)

    for dir in joinpath.(path_to_save_folder, ["data", "plots"])
        ispath(dir) || mkdir(dir)
    end
end

function ParsingData(path_to_file::String,
                      path_to_save::String)
    # Читаем исходник
    source_data = hton.(reinterpret(UInt32, read(path_to_file)))

    # Разбор 32-битных слов на компоненты
    control_bit = UInt8.((source_data .>> 31) .& 0x01)
    channel     = UInt8.((source_data .>> 25) .& 0x3F)
    delays      = UInt64.((source_data .>> 10) .& 0x7FFF)
    impulse_num = UInt64.(source_data .& 0x3FF)

    # Обработка переполнений impulse_num
    overflow_idxs = findall(==(1), control_bit)
    overflow_diff = Statistics.diff(overflow_idxs) .- 1
    overflow_vec  = collect(0 : length(overflow_diff)-1)
    overflow = Vector{UInt64}(undef, sum(overflow_diff))

    # АНАЛОГ repelem. МОДИФИЦИРОВАТЬ
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

    # Выбрасывание из массивов неинформативных данных через маску
    mask = ones(Bool, length(control_bit))
    mask[overflow_idxs] .= false

    channel     = channel[mask]
    delays      = delays[mask]
    impulse_num = impulse_num[mask]

    # Корректировка номеров импульсов и приведение задержек в корректную размерность
    impulse_num .+= overflow
    delays .*= 64

    # Формирование массива arrivals
    time_between = maximum(delays)
    arrivals = (delays .+ impulse_num .* UInt64(time_between)) .* 1e-12
    arrivals .-= arrivals[1]
    delays   = delays .* 1e-12

    # Сортировка чтобы исключить забегающие вперед фотоны
    perm    = sortperm(arrivals)
    delays  = delays[perm]
    channel = channel[perm]

    # Коррекция на длины плеч

    first_dec_idx  = findall(==(0), channel)
    second_dec_idx = findall(==(1), channel)

    first_dec_delay  = delays[first_dec_idx]
    second_dec_delay = delays[second_dec_idx]

    first_bins  = UInt32((maximum(first_dec_delay) - minimum(first_dec_delay)) * 1e12)
    second_bins = UInt32((maximum(second_dec_delay) - minimum(second_dec_delay)) * 1e12)

    decay_first  = fit(Histogram, first_dec_delay, nbins=first_bins)
    decay_second = fit(Histogram, second_dec_delay, nbins=second_bins)

    decay_first_max  = findmax(decay_first.weights)
    decay_second_max = findmax(decay_second.weights)

    decay_first_edges  = collect.(decay_first.edges)[1]
    decay_second_edges = collect.(decay_second.edges)[1]

    delays[first_dec_idx]  .-= decay_first_edges[decay_first_max[2]]
    delays[second_dec_idx] .-= decay_second_edges[decay_second_max[2]]

    mask = ones(Bool, length(delays))
    idx = findall(<(0), delays)
    mask[idx] .= false

    delays   = Float64.(delays[mask])
    arrivals = Float64.(arrivals[mask])
    channel  = UInt8.(channel[mask])

    parsed_data = ParsedData(channel, delays, arrivals)
    JLD2.save_object(path_to_save, parsed_data)
end

function MLP(params, N, T)
    tau_F, tau_D, A2, bkg = params

    if tau_F < 0 || tau_D < 0 || A2 < 0 || bkg < 0 || tau_F > tau_D
        return Inf
    end

    L = @. 1/tau_F * exp(-T/tau_F) + A2/tau_D * exp(-T/tau_D) + bkg
    L ./= sum(L)
    
    return -sum(N .* log.(L))
end