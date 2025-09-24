include("data_structure.jl")

using JLD2, Optim

function BinningData(path_to_parsed_data::String,
                    path_to_save::String,
                    bin_size::Float64)
    parsed_data = JLD2.load_object(path_to_parsed_data)

    arrivals = parsed_data.arrivals

    nums_of_bin = UInt32(arrivals[end] ÷ bin_size)
    nums_of_arrivals = length(arrivals)

    intensity = Vector{UInt32}(undef, nums_of_bin)
    index_min = Vector{UInt32}(undef, nums_of_bin)
    index_max = Vector{UInt32}(undef, nums_of_bin)

    count = 0
    for i in 1:nums_of_bin
        if i == 1
            index_min[i] = 1
            count = 1
        else
            index_min[i] = index_max[i-1] + 1
            count = index_min[i]
        end

        while arrivals[count] < i*bin_size
            count += 1
            if count + 1 > nums_of_arrivals
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

    binned_data = BinnedData(intensity, index_min, index_max)
    JLD2.save_object(path_to_save, binned_data)
end

function FLID_horizontal(
    path_to_parsed_data::String,
    path_to_save::String
)
    parsed_data = JLD2.load_object(path_to_parsed_data)

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
            p -> MLP(p, decay, decay_time),
            x0
        )
        tau_F[i], tau_D[i], A2[i], bkg[i] = Optim.minimizer(result)

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

    FLID_data = FLIDData(tau_F, tau_D, A2, bkg, I)
    JLD2.save_object(path_to_save, FLID_data)
end