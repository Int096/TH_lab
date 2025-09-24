using JLD2, CairoMakie, Optim, StatsBase, FHist, Colors

begin
    filename = "ProcessingResults/1_ms_per_bin/525/66/data/"

    binned_data = JLD2.load_object(
        filename*"binned_data.jld2"
    );
    parsed_data = JLD2.load_object(
        filename*"parsed_data.jld2"
    );

    intensity = binned_data.intensity;
    index_min = binned_data.index_min;
    index_max = binned_data.index_max;

    delays = parsed_data.delays * 1e9; 
end

begin
    nums_of_levels = 100

    intensity_borders = range(
        start = minimum(intensity),
        stop = maximum(intensity),
        length = nums_of_levels + 1
    )
end

nums_of_levels = 100
intensity_borders = range(
    start = minimum(intensity),
    stop = maximum(intensity),
    length = nums_of_levels+1
)

function MLP(params, N, T)
    tau_F, tau_D, A2, bkg = params

    if tau_F < 0 || tau_D < 0 || A2 < 0 || bkg < 0 || tau_F > tau_D
        return Inf
    end

    L = @. 1/tau_F * exp(-T/tau_F) + A2/tau_D * exp(-T/tau_D) + bkg
    L ./= sum(L)
    
    return -sum(N .* log.(L))
end

av_int_levels = (intensity_borders[1:end-1] + intensity_borders[2:end])/2
av_int_levels2 = intensity_borders[1:end-1]

array_cell = Vector{Vector{Int}}(undef, nums_of_levels)
for i in 1:nums_of_levels
    if i == 1
        kk = findall(
            intensity .<= intensity_borders[i+1]
        )
    else
        kk = findall(
            (intensity .<= intensity_borders[i+1])
            .&
            (intensity .> intensity_borders[i])
        )
    end
    array_cell[i] = kk
end

indexes_cell = Vector{Vector{Int}}(undef, nums_of_levels)
for i in 1:nums_of_levels
    println("i = $i")
    intensity_indexes = array_cell[i]
    #println("size(intensity_indexes) = $(size(intensity_indexes))")
    photon_indexes = []
    for j = intensity_indexes
        photon_indexes = vcat(
            photon_indexes,
            index_min[j] : index_max[j]
        )
    end
    indexes_cell[i] = photon_indexes
end

nums_of_ph_per_fit = 1000
count = 0
fit_initial_tau = range(0, 10, nums_of_levels)
fit_initial_int = (maximum(av_int_levels2)/10) .* fit_initial_tau


delays_data = []

for i in 1:nums_of_levels
    println("i = $i")

    indexes_in_level = indexes_cell[i]
    nm_ph = length(indexes_in_level)
    num_of_fits = nm_ph ÷ nums_of_ph_per_fit
    for j in 1:num_of_fits
        if j != num_of_fits
            data_indexes = indexes_in_level[
                ((j-1)*nums_of_ph_per_fit+1):j*nums_of_ph_per_fit
            ]
            data = delays[data_indexes]
            delays_data = data
        else 
            if j == 1
                data_indexes = indexes_in_level[
                    1:end
                ]
            else
                data_indexes = indexes_in_level[
                    ((j-1)*nums_of_ph_per_fit):end
                ]
            end
            data = delays[data_indexes]
        end

        hist_obj = fit(Histogram, data, 0.5:1:ceil(maximum(data)))
        decay = hist_obj.weights ./ sum(hist_obj.weights)
        decay_time = (hist_obj.edges[1:end-1] .+ hist_obj.edges[2:end])./2

        x0 = [av_int_levels2[i]/(maximum(av_int_levels2)*10) 1e-3 1e-5 100]


        result = Optim.optimize(
            p -> MLP(p, decay, decay_time),
            x0
        )
        println(Optim.minimizer(result))
    end
end

begin
    hist_obj = fit(Histogram, delays_data, 0.5:1:ceil(maximum(delays_data)))
        decay = hist_obj.weights ./ sum(hist_obj.weights)
        decay_time = (hist_obj.edges[1:end-1] .+ hist_obj.edges[2:end])./2

        fig = Figure()
        ax = Axis(fig[1,1])
        scatter!(decay_time, decay)

        save("FLID.png", figs)
end

