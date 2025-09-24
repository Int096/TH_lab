using JLD2, FHist, Optim, Statistics, CairoMakie, LaTeXStrings

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

folders = ["525", "525 PPD"]
path_to_save = "ProcessingResults/1_ms_per_bin/tau_F_distribution/"


function GenerateData(folders, path_to_save)
    for (i, folder) in enumerate(folders)
        path_to_data = "ProcessingResults/1_ms_per_bin/"*folder
        files = readdir(path_to_data)

        println("$i: $path_to_data")

        taus_mean = []
        taus_std = []

        sum = 0


        for (j, file) in enumerate(files)
            path_to_file = joinpath(path_to_data, file)
            println("\t$j: $path_to_file")

            GC.gc()

            parsed_data = JLD2.load_object(
                path_to_file*"/data/parsed_data.jld2"
            )
            binned_data = JLD2.load_object(
                path_to_file*"/data/binned_data.jld2"
            )

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
                sum = 0 
                #println("i = $i")
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
                    
                    #taus = params[1] .* ones(length(indexes_of_photons))
                    taus = params[1]
                    tau_F = vcat(tau_F, taus)

                    sum += taus * length(indexes_of_photons)
                end
            end

            taus_mean = vcat(taus_mean, sum / length(delays))
            taus_std  = vcat(taus_std, std(tau_F))
        end

        TauDistr = TausDistribution(
            folder,
            taus_mean,
            taus_std
        )
        JLD2.save_object(
            path_to_save*folder*"_old_pdf_data.jld2",
            TauDistr
        )
    end
end

function Plotter(
    folders,
    path_to_data,
    path_to_plot
)

    fig = Figure(
        size = (1200, 600),
        fontsize = 24
    )
    ax = Axis(
        fig[1, 1],
        xlabel = L"\text{Number of dots}",
        ylabel = L"Mean $\tau_F$, [ns]",
        title = "Распределение средних времен жизни; старый алгоритм (без нормировки)",
        limits = (0, 101, 0, 15)
        )

    colors = [:blue, :red]
    for (i, folder) in enumerate(folders) 
        println(path_to_data*"/"*folder*"_old_pdf_data.jld2")
        data = JLD2.load_object(path_to_data*"/"*folder*"_old_pdf_data.jld2")

        scatter!(
            ax,
            1:length(data.taus_mean),
            data.taus_mean,
            color = colors[i],
            label = folder*" (mean = $(round(mean(data.taus_mean);digits=2))±$(round(std(data.taus_mean); digits=2)))"
        )
    end
    axislegend()
    save(path_to_plot, fig)
end

GenerateData(folders, path_to_save)
Plotter(folders, path_to_save, 
        joinpath(path_to_save, "TAU_mean_distribution_old_pdf.png")
)