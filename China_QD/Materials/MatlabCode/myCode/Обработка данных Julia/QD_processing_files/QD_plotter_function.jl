using CairoMakie
using SpecialFunctions

function PlotterTrace(save_path::String, 
                      intensity::Vector{UInt32}, 
                      bin_size::Float64
                      )
    bin_size_ms = bin_size*1000

    fig = Figure(size=(1500, 500), fontsize=24)
    ax = Axis(fig[1, 1],
              xlabel="Time, [s]",
              ylabel="Intensity, [counts/$bin_size_ms ms]",
              title="55 QD 622 nm in PPD",
              limits=(0, length(intensity)*bin_size, 0, 1.1*maximum(intensity))
        )
    times = range(0.0; step=bin_size, length=length(intensity)) |> collect

    lines!(ax, times, intensity, color=:black)
    
    save("Trace.png", fig)
end

function PlotterCorrFunction(corr_data::CorrData)
    p = corr_data.g2 .- 1
    t = corr_data.times
    new_p = corr_data.new_p
    new_t = corr_data.new_t
    new_N_kor = corr_data.new_N_kor
    experiment_time = corr_data.experiment_time
    new_delta_t = corr_data.new_delta_t
    n_m = corr_data.n_m

    fig = Figure(resolution = (800, 600), backgroundcolor = :white)

    # Create axis with log scale on x-axis
    ax = Axis(fig[1, 1],
              xscale = log10,
              xlabel = "Time [s]",
              ylabel = "Autocorrelation function",
              xlabelsize = 24,
              ylabelsize = 24,
              xticklabelsize = 24,
              yticklabelsize = 24,
              xgridvisible = true,
              ygridvisible = true)

    lines!(ax, t, p)

    # Compute error bars (chi2inv equivalent in Julia using SpecialFunctions)
    chi2inv_lower = 0.5 * quantile(Chisq(2 * new_N_kor), 0.025)
    chi2inv_upper = 0.5 * quantile(Chisq(2 * new_N_kor .+ 2), 0.975)
    error_lower = (new_N_kor .- chi2inv_lower) ./ (new_delta_t * n_m^2 * experiment_time)
    error_upper = (chi2inv_upper .- new_N_kor) ./ (new_delta_t * n_m^2 * experiment_time)

    # Plot error bars
    errorbars!(ax, new_t, new_p, error_lower, error_upper, whiskerwidth = 6, color = :black)
    scatter!(ax, new_t, new_p, markersize = 6, color = :black, strokecolor = :black, strokewidth = 1.5)

    # Set axis limits
    ylims!(ax, low = 0.99 * minimum(p), high = maximum(p))
    xlims!(ax, low = 0.8 * minimum(t), high = 1.02 * maximum(t))

    ax.spinesvisible = true

    save("Corr.png", fig)
end