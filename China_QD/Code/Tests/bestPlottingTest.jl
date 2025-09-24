using JLD2, CairoMakie

begin
folder = "525"
ppd = "with PPD"
file = "100"


    bin_size = 0.001
    binned_data = JLD2.load_object(
        "ProcessingResults/1_ms_per_bin/525/1/data/binned_data.jld2"
    )
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
        color = :black
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