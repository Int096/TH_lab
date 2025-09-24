using JLD2, CairoMakie, StatsBase, Colors, Optim, FHist, Statistics

begin
    parsed_data = JLD2.load_object(
        "ProcessingResults/1_ms_per_bin/525 PPD/22/data/parsed_data.jld2"
    );
    binned_data = JLD2.load_object(
        "ProcessingResults/1_ms_per_bin/525 PPD/22/data/binned_data.jld2"
    );

    delays = parsed_data.delays * 1e9

    intensity = binned_data.intensity
    index_min = binned_data.index_min
    index_max = binned_data.index_max 
end

begin
    nums_of_levels = 100 

    intensity_border = range(
        start = minimum(intensity),
        stop = maximum(intensity),
        length = nums_of_levels+1
    )

    indexes_of_delays = Vector{Vector{Int}}(undef, nums_of_levels)
    delays_per_level = Vector{Vector{Int}}(undef, nums_of_levels)
    for i in 1 : nums_of_levels
        if i == 1
            index_of_intensity_with_borders = findall(
                intensity .<= intensity_border[i+1]
            )
        else 
            index_of_intensity_with_borders = findall(
                (intensity .<= intensity_border[i+1])
                .&
                (intensity .> intensity_border[i])
            )
        end
        indexes_of_delays[i] = index_of_intensity_with_borders

        photons_per_level = []
        for idx in index_of_intensity_with_borders
            photons_per_level = vcat(
                photons_per_level,
                index_min[idx] : index_max[idx]
            )
        end

        delays_per_level[i] = photons_per_level
    end
    delays_per_level
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

begin
    nums_of_photons_per_fit = 1000

    data = []

    tau_F = []
    tau_D = []
    A2    = []
    bkg   = []

    I = []

    av_ing_levels = (intensity_border[1:end-1] + intensity_border[2:end])/2
    
    for i in 1 : nums_of_levels
        indexes_in_level = delays_per_level[i]
        count = 0

        nums_of_ph = length(indexes_in_level)
        nums_of_fit = nums_of_ph ÷ nums_of_photons_per_fit

        for j = 1 : nums_of_fit 
            count += 1
            if j != nums_of_fit
                data_idx = indexes_in_level[
                    ((j-1)*nums_of_photons_per_fit+1):j*nums_of_photons_per_fit
                ]
                data = delays[data_idx]
            else
                if j == 1
                    data_indexes = indexes_in_level[
                        1:end
                    ]
                else
                    data_indexes = indexes_in_level[
                        ((j-1)*nums_of_photons_per_fit):end
                    ]
                end
                data = delays[data_indexes]
            end

            #hist_obj = fit(Histogram, 0.5 : (ceil(maximum(data))-0.5))
            #T = (collect.(hist_obj.edges)[1])[1:end-1] .+ 0.5
            #N = hist_obj.weights ./ sum(hist_obj.weights)

            hist_obj = Hist1D(data; binedges=0.5:(ceil(maximum(data))-0.5))
            T = bincenters(hist_obj)
            N = bincounts(hist_obj)

            x0 = [10., 100., 1., 1e-3]
            result = Optim.optimize(
                p -> MLP(p, N, T),
                x0
            )
            params = Optim.minimizer(
                result 
            )

            tau_F = vcat(tau_F, params[1])
            tau_D = vcat(tau_D, params[2])
            A2    = vcat(A2, params[3])
            bkg   = vcat(bkg, params[4])

            fn = @. 1/tau_F[end] * exp(-T/tau_F[end]) + A2[end]/tau_D[end] * exp(-T/tau_D[end]) + bkg[end]
            fn_F = sum(@. 1/tau_F[end] * exp(-T/tau_F[end])) / sum(fn)
            fn_D = sum(@. A2[end]/tau_D[end] * exp(-T/tau_D[end])) / sum(fn)

            I = vcat(I, av_ing_levels[i] * (fn_D + fn_F))
        end
    end
end


begin
matlab_parula = [
    RGB(0.2081, 0.1663, 0.5292),
    RGB(0.2116, 0.1898, 0.5777),
    RGB(0.2123, 0.2138, 0.6270),
    RGB(0.2081, 0.2386, 0.6771),
    RGB(0.1959, 0.2645, 0.7279),
    RGB(0.1707, 0.2919, 0.7792),
    RGB(0.1253, 0.3242, 0.8303),
    RGB(0.0591, 0.3598, 0.8683),
    RGB(0.0117, 0.3875, 0.8820),
    RGB(0.0060, 0.4086, 0.8828),
    RGB(0.0165, 0.4266, 0.8786),
    RGB(0.0329, 0.4430, 0.8720),
    RGB(0.0498, 0.4586, 0.8641),
    RGB(0.0629, 0.4737, 0.8554),
    RGB(0.0723, 0.4887, 0.8467),
    RGB(0.0779, 0.5040, 0.8384),
    RGB(0.0793, 0.5200, 0.8312),
    RGB(0.0749, 0.5375, 0.8263),
    RGB(0.0641, 0.5570, 0.8240),
    RGB(0.0488, 0.5772, 0.8228),
    RGB(0.0343, 0.5966, 0.8199),
    RGB(0.0265, 0.6137, 0.8135),
    RGB(0.0239, 0.6287, 0.8038),
    RGB(0.0231, 0.6418, 0.7913),
    RGB(0.0228, 0.6535, 0.7768),
    RGB(0.0267, 0.6642, 0.7607),
    RGB(0.0384, 0.6743, 0.7436),
    RGB(0.0590, 0.6838, 0.7254),
    RGB(0.0843, 0.6928, 0.7062),
    RGB(0.1133, 0.7015, 0.6859),
    RGB(0.1453, 0.7098, 0.6646),
    RGB(0.1801, 0.7177, 0.6424),
    RGB(0.2178, 0.7250, 0.6193),
    RGB(0.2586, 0.7317, 0.5954),
    RGB(0.3022, 0.7376, 0.5712),
    RGB(0.3482, 0.7424, 0.5473),
    RGB(0.3953, 0.7459, 0.5244),
    RGB(0.4420, 0.7481, 0.5033),
    RGB(0.4871, 0.7491, 0.4840),
    RGB(0.5300, 0.7491, 0.4665),
    RGB(0.5709, 0.7485, 0.4504),
    RGB(0.6099, 0.7473, 0.4356),
    RGB(0.6473, 0.7456, 0.4218),
    RGB(0.6834, 0.7435, 0.4088),
    RGB(0.7184, 0.7411, 0.3965),
    RGB(0.7525, 0.7384, 0.3847),
    RGB(0.7858, 0.7356, 0.3735),
    RGB(0.8185, 0.7327, 0.3629),
    RGB(0.8507, 0.7298, 0.3529),
    RGB(0.8824, 0.7269, 0.3435),
    RGB(0.9139, 0.7241, 0.3347),
    RGB(0.9450, 0.7215, 0.3265),
    RGB(0.9739, 0.7193, 0.3190),
    RGB(0.9932, 0.7176, 0.3127),
    RGB(0.9990, 0.7164, 0.3075),
    RGB(0.9955, 0.7157, 0.3035),
    RGB(0.9880, 0.7154, 0.3006),
    RGB(0.9789, 0.7154, 0.2989),
    RGB(0.9697, 0.7157, 0.2983),
    RGB(0.9626, 0.7162, 0.2988),
    RGB(0.9589, 0.7169, 0.3003),
    RGB(0.9598, 0.7177, 0.3028),
    RGB(0.9661, 0.7186, 0.3062),
    RGB(0.9763, 0.7195, 0.3104),
    RGB(0.9879, 0.7205, 0.3152),
    RGB(0.9982, 0.7215, 0.3205),
    RGB(1.0000, 0.7225, 0.3262)
]

folder = "525"
ppd = ""
file = "12"
    h2d = Hist2D(
        (tau_F, I),
        nbins=(150, 150),
        overflow = true
        )

    fig, ax, p = plot(
        h2d,
        figure = (
            size=(800, 600),
            fontsize = 24,
            font="Times New Roman",
            ),
        axis = (
            backgroundcolor = RGB(0., 0., 0.49),
            xlabel = "Fast time, [ns]",
            ylabel = "Yield",
        ),
        colormap = matlab_parula,
    )
    cbar = Colorbar(fig[1, 2], p, width=15, tickwidth=2, tickalign=1.5)
    
    # Добавляем вертикальную подпись слева от колорбара
    Label(fig[1, 2, Left()], "Occurence",
        rotation=π/2,  # Поворот на 90 градусов
        color=:black  # Отступы (left, right, bottom, top)
    )
    text!(
        ax,
        "QD №$file\n$(split(folder)[1]) nm $ppd",
        position = (0.95, 0.95), 
        space = :relative,
        align = (:right, :top),
        color = :white,
        fontsize = 20
    )
    save("FLID.png", fig)
    fig
end

begin
    hist_obj = Hist1D(tau_F; binedges=0.5:(ceil(maximum(tau_F))-0.5))
    T = bincenters(hist_obj)
    N = bincounts(hist_obj)

    println(mean(tau_F))
    println(std(tau_F))

    scatter(T, N)
end