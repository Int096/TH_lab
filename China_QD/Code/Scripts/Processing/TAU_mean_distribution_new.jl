using JLD2, CairoMakie, LaTeXStrings, Statistics, FHist, Optim

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
path_to_data = joinpath("ProcessingResults", "1_ms_per_bin", "tau_F_distribution")
path_to_plot = joinpath(path_to_data, "TAU__mean_distribution_new.png")

function GenerateData()
    data_folders = ["525 PPD", "525"]

    path_to_Result = joinpath("ProcessingResults", "1_ms_per_bin", )

    for (i, folder) in enumerate(data_folders) 
        path_to_data_folder = joinpath(path_to_Result, folder)
        files = readdir(path_to_data_folder)

        println("$i: $path_to_data_folder")

        taus_mean = []
        taus_std  = []

        for (j, file) in enumerate(files) 
            path_to_file = joinpath(path_to_data_folder, file)
            println("\t$j: $path_to_file")

            parsed_data = JLD2.load_object(
                path_to_file*"/data/parsed_data.jld2"
            )

            delays = parsed_data.delays * 1e9
            arrivals = parsed_data.arrivals

            num_ph_per_bin = 1000
            N = length(delays) ÷ num_ph_per_bin

            tau_F = Vector{Float64}(undef, N)

            for i in 1:N
                k = ((i-1)*num_ph_per_bin+1) : i*num_ph_per_bin

                data = delays[k]

                hist_obj = Hist1D(data; binedges=0.5:(ceil(maximum(data))-0.5))
                T = bincenters(hist_obj)
                N = bincounts(hist_obj) ./ sum(bincounts(hist_obj))

                x0 = [10., 100., 1., 1e-3]
                result = Optim.optimize(
                    p -> MLP(p, N, T),
                    x0
                )
                tau_F[i], tau_D, A2, bkg = Optim.minimizer(result)
            end

            taus_mean = vcat(taus_mean, mean(tau_F))
            taus_std  = vcat(taus_std, std(tau_F))
        end
        TauDistr = TausDistribution(
            folder,
            taus_mean,
            taus_std
        )

        JLD2.save_object(
            "ProcessingResults/1_ms_per_bin/tau_F_distribution/"*folder*"_data.jld2",
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
        title = "Распределение средних времен жизни; новый алгоритм (без нормировки)",
        limits = (0, 101, 0, 15)
        )

    colors = [:blue, :red]
    for (i, folder) in enumerate(folders) 
        println(path_to_data*"/"*folder*"_data.jld2")
        data = JLD2.load_object(path_to_data*"/"*folder*"_data.jld2")

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

#GenerateData()
Plotter(folders, path_to_data, path_to_plot)
