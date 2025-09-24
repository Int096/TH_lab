using JLD2, Optim, StatsBase, CairoMakie, QuadGK

function MLP(params, int_centers, occurence)
    S(Y, alpha, beta) = @. beta^(1/alpha) * (1 - Y)^(1/alpha) * Y^(-1/alpha)
    p(Y, alpha, beta, S_bar, sigma) = @. beta^(1/alpha) / alpha * (1 - Y)^((1-alpha)/alpha) * Y^(-(1+alpha)/alpha) * exp(-1/(2*sigma^2) * (S(Y, alpha, beta) - S_bar)^2) 
    p_norm(Y, alpha, beta, S_bar, sigma) = @. p(Y, alpha, beta, S_bar, sigma) / quadgk(Y -> p(Y, alpha, beta, S_bar, sigma), 0, 1)[1]

    alpha = 10.0
    beta, S_bar, sigma = params

    if beta < 0 || S_bar < 0 || sigma < 0 || alpha < 0
        return Inf
    end

    n_mean = mean(int_centers)
    Y_mean = quadgk(Y -> Y*p_norm(Y, alpha, beta, S_bar, sigma), 0, 1)[1]
    N = n_mean / Y_mean

    L = @. beta^(1/alpha) * N^int_centers / factorial(big(int_centers)) * quadgk(Y -> (1-Y)^((1-alpha)/alpha) * Y^(((int_centers-1)*alpha-1)/alpha) * exp(-1/(2*sigma^2) * (S(Y, alpha, beta) - S_bar)^2), 0, 1)[1]
    L ./= sum(L)

    return -sum(occurence .* log.(L))
end

function Q(params, int_centers, occurence)
    S(Y, alpha, beta) = beta^(1/alpha) * (1 - Y)^(1/alpha) * Y^(-1/alpha)
    p(Y, alpha, beta, S_bar, sigma) = beta^(1/alpha) / alpha * (1 - Y)^((1-alpha)/alpha) * Y^(-(1+alpha)/alpha) * exp(-1/(2*sigma^2) * (S(Y, alpha, beta) - S_bar)^2) 
    p_norm(Y, alpha, beta, S_bar, sigma) = p(Y, alpha, beta, S_bar, sigma) / quadgk(Y -> p(Y, alpha, beta, S_bar, sigma), 0, 1)[1]

    alpha = 10.0
    beta, S_bar, sigma = params

    n_mean = mean(int_centers)
    Y_mean = quadgk(Y -> Y*p_norm(Y, alpha, beta, S_bar, sigma), 0, 1)[1]
    N = n_mean / Y_mean

    Q = @. beta^(1/alpha) * N^int_centers / factorial(big(int_centers)) * quadgk(Y -> (1-Y)^((1-alpha)/alpha) * Y^(((int_centers-1)*alpha-1)/alpha) * exp(-1/(2*sigma^2) * (S(Y, alpha, beta) - S_bar)^2), 0, 1)[1]
    return Q
end

folder = "525"
file = "89"

binned_data = JLD2.load_object("../results/"*folder*"/"*file*"/data/binned_data.jld2")
intensity = binned_data.intensity

min, max = minimum(intensity), maximum(intensity)

hist_obj = fit(Histogram, intensity, min-0.5:max+0.5)

int_centers = min : max
occurence = hist_obj.weights ./ sum(hist_obj.weights)

##########
#
#
##########

x0 = [100., 10., 10.]
options = Optim.Options(
                        iterations = 10_000,
                        g_tol = 1e-12
                       )
opt = Optim.Options(g_tol = 1e-12,
                             iterations = 10,
                             store_trace = true,
                             show_trace = false,
                             show_warnings = true)
result = Optim.optimize(
                        p -> MLP(p, int_centers, occurence),
                        x0,
                        opt
                       )
params = Optim.minimizer(result)

println(result)
println(params)

##########

fig = Figure(
             size = (1200, 600),
             fontsize = 24
            )
ax = Axis(
          fig[1, 1],
          xlabel = "Intensity",
          ylabel = "Occurence",
          limits = (0, maximum(int_centers), 0, 1.05*maximum(occurence))
         )
ax1 = Axis(
          fig[2, 1],
          xlabel = "Intensity",
          ylabel = "Occurence" 
         )

barplot!(
         ax,
         int_centers, 
         occurence,
         color = :gray,
         gap = 0,
         alpha = 1.0,
         strokewidth = 0
        )
scatter!(
         ax,
         int_centers,
         occurence,
         color = :black,
         label = "Experiment"
        )
lines!(
       ax1,
       int_centers,
       Q(params, int_centers, occurence),
       color = :red,
       label = "Best fit: beta = $(round(params[1]; digits=2)), S_bar = $(round(params[2];digits=2)), sigma = $(round(params[3]; digits=2))"
      )

axislegend()

save("Occurence_and_Intensity.png", fig)
fig
