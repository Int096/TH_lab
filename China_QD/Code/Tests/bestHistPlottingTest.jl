using JLD2, CairoMakie, StatsBase

begin
folder = "525"
ppd = "with PPD"
file = "100"


    bin_size = 0.001
    binned_data = JLD2.load_object(
        "ProcessingResults/1_ms_per_bin/525/66/data/binned_data.jld2"
    )
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
        ylabel="Occurence, normalized",
        limits=(0, maximum(intensity), 0, 1.01*maximum(weights))
        )
    barplot!(
        ax,
        centers,
        weights,
        color=:black,
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
        fig
end