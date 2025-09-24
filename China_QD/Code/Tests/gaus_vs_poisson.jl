using CairoMakie

function Gauss(x, x_mean, sigma)
    return exp(-1/(2*sigma^2) * (x - x_mean)^2)
end

let x =  0:0.001:5,
    x_discr = 0:1:5

    gauss = Gauss.(x, 2, 1)

    fig = Figure(
        size = (800, 600), 
        fontsize = 24,
        lims = (0, 5, 0, maximum(gauss))
    )
    ax = Axis(
        fig[1, 1],
        xlabel = "x", 
        ylabel = "y",
        title = "Gauss vs Poisson"
    )

    lines!(
        ax,
        x,
        gauss,
        color = :blue,
        label = "Gauss; mean = 2, std = 1"
    )

    axislegend()
    save("trash/gauss_vs_poisson.png", fig)
    fig

end