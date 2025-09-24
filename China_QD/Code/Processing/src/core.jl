using Statistics, JLD2, StatsBase, Optim

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

struct newFLIDStruct
    tau_F::Vector{Float64}
    I::Vector{Float64}
end

function MLP_biex(params, N, T)
    tau_F, tau_D, A2, bkg = params

    if tau_F < 0 || tau_D < 0 || A2 < 0 || bkg < 0 || tau_F > tau_D
        return Inf
    end

    L = @. 1/tau_F * exp(-T/tau_F) + A2/tau_D * exp(-T/tau_D) + bkg
    L ./= sum(L)
    
    return -sum(N .* log.(L))
end

function ParseDataFromBinToJLD(
    path_to_file_with_raw_data::String,
    path_to_file_with_parsed_data::String
    )
    # Чтение исходника из бинарного файла  

    source_data = hton.(reinterpret(UInt32, read(path_to_file_with_raw_data)))

    # Разбор 32-битных слов на компоненты 
    control_bit = UInt8.((source_data .>> 31) .& 0x01)
    channel     = UInt8.((source_data .>> 25) .& 0x3F)
    delays      = UInt64.((source_data .>> 10) .& 0x7FFF)
    impulse_num = UInt64.(source_data .& 0x3FF)

    # Обработка переполнений impulse_num
    overflow_idxs = findall(==(1), control_bit)
    overflow_diff = Statistics.diff(overflow_idxs) .- 1
    overflow_vec  = collect(0 : length(overflow_diff)-1)
    
    overflow = [v for (v, n) in zip(overflow_vec, overflow_diff) for _ in 1:n] .* 1023
    
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

    # Сохраняем данные
    parsed_data = ParsedData(
        channel,
        delays,
        arrivals
    )
    JLD2.save_object(
        path_to_file_with_parsed_data,
        parsed_data
    )
end

function BinningData(
    path_to_file_with_parsed_data::String,
    path_to_file_with_binned_data::String,
    bin_size::Float64
    )
    parsed_data = JLD2.load_object(
        path_to_file_with_parsed_data
    )

    arrivals = parsed_data.arrivals

    num_of_bin = UInt32(arrivals[end] ÷ bin_size)
    num_of_arrivals = length(arrivals)

    intensity = Vector{UInt32}(undef, num_of_bin)
    index_min = Vector{UInt32}(undef, num_of_bin)
    index_max = Vector{UInt32}(undef, num_of_bin)

    count = 0
    for i in 1:num_of_bin
        if i == 1
            index_min[i] = 1
            count = 1
        else
            index_min[i] = index_max[i-1] + 1
            count = index_min[i]
        end

        while arrivals[count] < i*bin_size
            count += 1
            if count + 1 > num_of_arrivals
                break
            end
        end

        if count == index_min[i]
            intensity[i] = 0
            index_max[i] = index_min[i] - 1
        else
            index_max[i] = count - 1
            intensity[i] = index_max[i] - index_min[i] + 1
        end
    end

    binned_data = BinnedData(
        intensity,
        index_min, 
        index_max
    )
    JLD2.save_object(
        path_to_file_with_binned_data,
        binned_data
    )
end

function newFLIDData(
    path_to_file_with_parsed_data::String, 
    path_to_file_with_newFLID_data::String
    )
    parsed_data = JLD2.load_object(
        path_to_file_with_parsed_data
    )

    delays = parsed_data.delays * 1e9
    arrivals = parsed_data.arrivals

    num_ph_per_bin = 1000
    N = length(delays) ÷ num_ph_per_bin

    tau_F = Vector{Float64}(undef, N)
    tau_D = Vector{Float64}(undef, N)
    A2    = Vector{Float64}(undef, N)
    bkg   = Vector{Float64}(undef, N)

    Rtime = Vector{Float64}(undef, N)
    I     = Vector{Float64}(undef, N)

    for i in 1:N 
        k = ((i-1)*num_ph_per_bin+1) : i*num_ph_per_bin

        Rtime[i] = arrivals[k[end]]

        data = delays[k]
        
        hist_obj = fit(Histogram, data, 0.5 : (ceil(maximum(data))-0.5))
        decay_time = (collect.(hist_obj.edges)[1])[1:end-1] .+ 0.5
        decay = hist_obj.weights ./ sum(hist_obj.weights)
        
        x0 = [10., 100., 1., 1e-3]
        result = Optim.optimize(
            p -> MLP_biex(p, decay, decay_time),
            x0,
            NelderMead(),
            Optim.Options(
                iterations = 10000
            )
        )
        tau_F[i], tau_D[i], A2[i], bkg[i] = Optim.minimizer(result)

        if !Optim.converged(result)
            println("ВСЕ СЛОМАЛОСЬ")
        end

        fn = @. 1/tau_F[i] * exp(-decay_time/tau_F[i]) + A2[i]/tau_D[i] * exp(-decay_time/tau_D[i]) + bkg[i]
        fn_F = sum(@. 1/tau_F[i] * exp(-decay_time/tau_F[i])) / sum(fn)
        fn_D = sum(@. A2[i]/tau_D[i] * exp(-decay_time/tau_D[i])) / sum(fn)

        if i > 1
            inten = num_ph_per_bin / (Rtime[i]-Rtime[i-1]) / 1000
        else
            inten = num_ph_per_bin / Rtime[i] / 1000
        end

        I[i] = inten * (fn_D + fn_F)
    end

    FLID_data = newFLIDStruct(
        tau_F, 
        I
    )

    JLD2.save_object(
        path_to_file_with_newFLID_data,
        FLID_data
    )

end