using JLD2, Optim, FHist, CairoMakie, LaTeXStrings, Distributions, Statistics

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

struct WeightedTau
    tau::Float64
    weight::UInt32
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

function MLP_gauss(params, hist_X, hist_Y)
    sigma_1, mu_1, A2, sigma_2, mu_2 = params

    if sigma_1 < 0 || mu_1 < 0 || A2 < 0 || sigma_2 < 0 || mu_2 < 0
        return Inf 
    end

    L = @. 1/sqrt(2*π*sigma_1^2) * exp(-(hist_X - mu_1)^2/(2*sigma_1)) + A2/sqrt(2*π*sigma_2^2) * exp(-(hist_X - mu_2)^2/(2*sigma_2))
    L ./= sum(L)

    return -sum(hist_Y .* log.(L))
end



function ProcessingData(folders, path_to_save)
    for (i, folder) in enumerate(folders)
        path_to_data = joinpath("ProcessingResults", "1_ms_per_bin", folder)
        files = readdir(path_to_data)

        println("$i/$(length(folders)): $path_to_data")

        dict_of_data = Dict{String, Vector{WeightedTau}}()

        for (j, file) in enumerate(files)
            path_to_file = joinpath(path_to_data, file)
            println("\tФайл $j/$(length(files)): $file")

            ###############################################
            #   Непосредственно обработка
            ###############################################

            parsed_data = JLD2.load_object(
                path_to_file*"/data/parsed_data.jld2"
            )
            binned_data = JLD2.load_object(
                path_to_file*"/data/binned_data.jld2"
            )

            # Распаковываем данные 
            delays = parsed_data.delays * 1e9

            intensity = binned_data.intensity 
            index_max = binned_data.index_max
            index_min = binned_data.index_min

            num_of_levels = 100 
            num_of_photons_per_fit = 1000

            intensity_border = range(
                start  = minimum(intensity),
                stop   = maximum(intensity),
                length = num_of_levels + 1
            )

            weighted_taus = Vector{WeightedTau}()

            # Итерируемся по уровням
            for m in 1:num_of_levels 
                if m == 1 
                    intensitys_in_borders = findall(
                        intensity .<= intensity_border[m+1]
                    )
                else 
                    intensitys_in_borders = findall(
                        (intensity .<= intensity_border[m+1])
                        .&
                        (intensity .> intensity_border[m])
                    )
                end

                indexes_of_photons = []
                for idx in intensitys_in_borders
                    indexes_of_photons = vcat(
                        indexes_of_photons,
                        index_min[idx] : index_max[idx]
                    )
                end

                num_of_fit = length(indexes_of_photons) ÷ num_of_photons_per_fit

                for k = 1:num_of_fit
                    if k != num_of_fit
                        data_idx = indexes_of_photons[
                            ((k-1)*num_of_photons_per_fit+1):k*num_of_photons_per_fit
                        ]
                        data = delays[data_idx]
                    else
                        if k == 1
                            data_idx = indexes_of_photons[
                                1:end
                            ]
                        else
                            data_idx = indexes_of_photons[
                                ((k-1)*num_of_photons_per_fit):end
                            ]
                        end
                        data = delays[data_idx]
                    end

                    hist_obj = Hist1D(
                        data; 
                        binedges = 0.5 : (ceil(maximum(data))-0.5)
                    ) |> normalize
                    hist_X = bincenters(hist_obj)
                    hist_Y = bincounts(hist_obj)

                    fit_result = Optim.optimize(
                        p -> MLP(p, hist_Y, hist_X),
                        [10., 100., 1., 1e-3]
                    )
                    fit_parametrs = Optim.minimizer(fit_result)

                    if !Optim.converged(fit_result)
                        println("AAAAAAAAAAAAAAAAAAAAAAA")
                    end

                    taus = WeightedTau(fit_parametrs[1], length(indexes_of_photons))

                    push!(weighted_taus, taus)
                end
            end

            dict_of_data[file] = weighted_taus
        end

        save(path_to_save*"/"*folder*"_old_pdf_data.jld2", "dict_of_data", dict_of_data)
    end
end

function Plotter(folders, path_to_save)
    fig = Figure(
        size = (2000, 1000),
        fontsize = 34
    )
    gl = GridLayout(fig[1, 1])
    Label(
        gl[0, 1:6],
        "Распределение средних времен жизни\nВзвешенное среднее по количеству фотонов в уровне",
        fontsize = 28,
        font = :bold,
        halign = :center
    )
    ax = Axis(
        gl[1, 1:4],
        xlabel = L"\text{Number of dots}",
        ylabel = L"Mean $\tau_F$, [ns]",
        #title = "Распределение средних времен жизни;\nВзвешенное среднее по количеству фотонов в уровне",
        limits = (0, 120, 0, 15)
    )
    ax_hist_525 = Axis(
        gl[1, 5],
        ylabelvisible = false,     
        yticklabelsvisible = false,  
        yticksvisible = false,
        xlabelvisible = false,     
        xticklabelsvisible = false,  
        xticksvisible = false   
    )
    ax_hist_525_PPD = Axis(
        gl[1, 6],
        ylabelvisible = false,     
        yticklabelsvisible = false,  
        yticksvisible = false,
        xlabelvisible = false,     
        xticklabelsvisible = false,  
        xticksvisible = false   
    )
    colgap!(gl, 0)
    rowgap!(gl, 20)

    Y_max = 0

    colors = [:blue, :red]
    for (i, folder) in enumerate(folders)
        println("Отрисовка данных из каталога $folder ($i/$(length(folders)))")

        path_to_data = joinpath("ProcessingResults", "1_ms_per_bin", folder)
        files = readdir(path_to_data)

        loaded_data = load(path_to_save*"/"*folder*"_old_pdf_data.jld2")
        data = loaded_data["dict_of_data"]

        tau_means = Vector{Float64}()
        number_of_dots = Vector{Int32}()

        for (j, file) in enumerate(files)
            sum = 0
            W = 0
            for taus in data[file]
                sum += taus.tau * taus.weight
                W += taus.weight
            end

            push!(tau_means, sum/W)
            push!(number_of_dots, parse(Int32, file))
        end

        scatter!(
            ax,
            number_of_dots,
            tau_means,
            color = colors[i],
            markersize = 18,
            label = folder*" (μ = $(round(mean(tau_means);digits=2)), σ = $(round(std(tau_means);digits=2)))"
        )

        hist_obj = Hist1D(
            tau_means;
            binedges = 0:15
        ) |> normalize
        hist_Y = bincounts(hist_obj)
        hist_X = bincenters(hist_obj)

        Y_max = max(Y_max, maximum(hist_Y))
    
        if folder == "525"
            barplot!(
                ax_hist_525,
                hist_X,
                hist_Y,
                color = colors[i],
                gap = 0,
                direction = :x
            )
        else
            barplot!(
                ax_hist_525_PPD,
                hist_X,
                hist_Y,
                color = colors[i],
                gap = 0,
                direction = :x
            )
        end
    end

    xlims!(ax_hist_525, 0, 1.08*Y_max)
    ylims!(ax_hist_525, 0, 15)
    xlims!(ax_hist_525_PPD, 0, 1.08*Y_max)
    ylims!(ax_hist_525_PPD, 0, 15)

    axislegend(ax, position = (:right, :bottom))
    save(path_to_save*"/plot.png", fig)
    fig
end

function PlotterHist(folders, path_to_save)
    for (i, folder) in enumerate(folders)
        loaded_data = load(path_to_save*"/"*folder*"_old_pdf_data.jld2")
        data = loaded_data["dict_of_data"]

        path_to_data = joinpath("ProcessingResults", "1_ms_per_bin", folder)
        files = readdir(path_to_data)

        for (j, file) in enumerate(files)
            tau = [weighted_tau.tau for weighted_tau in data[file]]
            weight = [weighted_tau.weight for weighted_tau in data[file]]

            hist_obj = Hist1D(
                tau;
                weights = weight,
                binedges = 0 : ceil(maximum(tau))
            ) |> normalize
            hist_X = bincenters(hist_obj)
            hist_Y = bincounts(hist_obj)

            fit_result = Optim.optimize(
                p -> MLP_gauss(p, hist_X, hist_Y),
                [0.01, 2., 1e-3, 0.1, 11],
                Optim.Options(
                    g_tol = 1e-12
                )
            )
            fit_params = Optim.minimizer(fit_result)

            println(fit_result)

            sigma_1, mu_1, A2, sigma_2, mu_2 = fit_params
            fn = @. 1/sqrt(2*π*sigma_1^2) * exp(-(hist_X - mu_1)/(2*sigma_1)) + A2/sqrt(2*π*sigma_2^2) * exp(-(hist_X - mu_2)/(2*sigma_2))

            fig = Figure(
                size = (800, 600),
                fontsize = 24
            )
            ax = Axis(
                fig[1,1],
                xlabel = "Fast time, [ns]",
                ylabel = "Occurence",
                limits = (0, maximum(hist_X), 0, 1.08*maximum(hist_Y))
            )
            barplot!(
                ax, 
                hist_X,
                hist_Y,
                color = :gray,
                gap = 0
            )
            scatter!(
                ax, 
                hist_X,
                hist_Y, 
                color = :black
            )
            lines!(
                ax, 
                hist_X,
                fn, 
                color = :red
            )
           
            save(path_to_save*"/"*folder*"_"*file*".png", fig)
        end
    end
end

#ProcessingData(
#    ["525", "525 PPD"],
#    "ProcessingResults/1_ms_per_bin/tau_F_distribution"
#)

#Plotter(
#    ["525", "525 PPD"],
#    "ProcessingResults/1_ms_per_bin/tau_F_distribution"
#)

PlotterHist(
    ["525"],
    "ProcessingResults/1_ms_per_bin/tau_F_distribution"
)