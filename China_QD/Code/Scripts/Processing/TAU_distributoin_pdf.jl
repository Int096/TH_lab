using JLD2, FHist, Optim, CairoMakie, LaTeXStrings, Statistics

struct ParsedData
    channel::Vector{UInt32}
    delays::Vector{Float64}
    arrivals::Vector{Float64}
end

struct TausDistribution 
    folder::String 
    taus_mean::Vector{Float64}
    taus_std::Vector{Float64}
end

struct BinnedData
    intensity::Vector{UInt32}
    index_min::Vector{UInt32}
    index_max::Vector{UInt32}
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

let folder = "525", file = "22"
    path_to_data = "ProcessingResults/1_ms_per_bin/"*folder*"/"*file*"/data/"

    parsed_data = JLD2.load_object(
        path_to_data*"parsed_data.jld2"
    );
    binned_data = JLD2.load_object(
        path_to_data*"binned_data.jld2"
    );

    delays = parsed_data.delays * 1e9

    intensity = binned_data.intensity 
    index_min = binned_data.index_min
    index_max = binned_data.index_max

    nums_of_level = 100
    nums_of_photons_per_fit = 1000

    intensity_border = range(
        start = minimum(intensity),
        stop = maximum(intensity),
        length = nums_of_level + 1
    )

    tau_F = []

    for i in 1:nums_of_level
        println(i)
        # Бины, попавшие в границы
        if i == 1
            intensity_in_borders = findall(
                intensity .<= intensity_border[i+1]
            )
        else
            intensity_in_borders = findall(
                (intensity .<= intensity_border[i+1])
                .&
                (intensity .> intensity_border[i])
            )
        end

        # Индексы фотонов
        indexes_of_photons = []
        for idx in intensity_in_borders
            indexes_of_photons = vcat(
                indexes_of_photons,
                index_min[idx] : index_max[idx]
            )
        end

        # Движемся по фотонам и фитим
        nums_of_fit = length(indexes_of_photons) ÷ nums_of_photons_per_fit

        for j = i:nums_of_fit 
            if j != nums_of_fit
                data_idx = indexes_of_photons[
                    ((j-1)*nums_of_photons_per_fit+1):j*nums_of_photons_per_fit
                ]
                data = delays[data_idx]
            else
                if j == 1
                    data_idx = indexes_of_photons[
                        1:end
                    ]
                else
                    data_idx = indexes_of_photons[
                        ((j-1)*nums_of_photons_per_fit):end
                    ]
                end
                data = delays[data_idx]
            end

            hist_obj = Hist1D(
                data;
                binedges = 0.5 : (ceil(maximum(data))-0.5)
            ) |> normalize
            T = bincenters(hist_obj)
            N = bincounts(hist_obj)

            result = Optim.optimize(
                p -> MLP(p, N, T),
                [10., 100., 1., 1e-3]
            )
            params = Optim.minimizer(
                result
            )
            
            taus = params[1] .* ones(length(indexes_of_photons))
            tau_F = vcat(tau_F, taus)

        end
    end


    hist_obj = Hist1D(
        tau_F;
        binedges = 0.5 :0.5: (ceil(maximum(tau_F))-0.5)
    ) |> normalize
    T = bincenters(hist_obj)
    N = bincounts(hist_obj)

    fig = Figure(
        size = (1000, 800),
        fontsize = 24
    )
    ax = Axis(
        fig[1, 1],
        xlabel = "Fast time, [ns]",
        ylabel = "Occurence",
        title = "Распределение времени жизни",
        limits = (0, maximum(T)+0.5, 0, 1.01*maximum(N))
    )

    barplot!(
        ax,
        T,
        N,
        gap = 0,
        color = :gray``
    )
    scatter!(
        ax,
        T,
        N,
        color = :black
    )

    fig
end