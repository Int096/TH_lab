include("data_structure.jl")

using JLD2, CairoMakie, Colors, FHist

function PlotterTrace(path_to_binned_data::String,
                      path_to_save::String,
                      bin_size::Float64,
                      folder::String,
                      file::String,
                      ppd::String)
    binned_data = JLD2.load_object(path_to_binned_data)
    intensity = binned_data.intensity

    fig = Figure(
            size = (1200, 600),
            fontsize = 24,
            font = "Times New Roman"
    )
    ax = Axis(
            fig[1, 1],
            xlabel = "Time, [s]",
            ylabel = "Intensity, [count/$(UInt32(bin_size*1000)) ms]",
            limits=(0, length(intensity)*bin_size, 0, 1.1*maximum(intensity))
    )

    times = range(0.0; step=bin_size, length=length(intensity)) |> collect

    lines!(
        ax,
        times,
        intensity,
        color = :gray
    )

    text!(
        ax,
        "QD №$file\n$(split(folder)[1]) nm $ppd",
        position = (0.95, 0.95), 
        space = :relative,
        align = (:right, :top)
    )

    save(path_to_save, fig)
end

function PlotterHist(path_to_binned_data::String,
                     path_to_save::String,
                     folder::String,
                     file::String,
                     ppd::String)

    binned_data = JLD2.load_object(path_to_binned_data)
    intensity = binned_data.intensity

    min_val, max_val = minimum(intensity), maximum(intensity)
    edges = min_val - 0.5 : 1.0 : max_val + 0.5
    hist_obj = fit(Histogram, intensity, edges)

    centers = min_val : 1 : max_val
    weights = hist_obj.weights ./ sum(hist_obj.weights)

    fig = Figure(
        size=(800, 600),
        fontsize=24
        )
    ax = Axis(
        fig[1, 1],
        xlabel="Intensity",
        ylabel="Occurence",
        limits=(0, maximum(intensity), 0, 1.01*maximum(weights))
        )
    barplot!(
        ax,
        centers,
        weights,
        color=:gray,
        strokewidth=0,
        gap=0,
        alpha=1.0
        )
    
    text!(
        ax,
        "QD №$file\n$(split(folder)[1]) nm $ppd",
        position = (0.95, 0.95), 
        space = :relative,
        align = (:right, :top)
    )
    save(path_to_save, fig)
end

function PlotterFLID(
    path_to_FLID_data::String,
    path_to_save::String,
    folder::String,
    file::String,
    ppd::String
)
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

    FLID_data = JLD2.load_object(path_to_FLID_data)
    tau_F = FLID_data.tau_F
    I = FLID_data.I

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
    save(path_to_save, fig)
end