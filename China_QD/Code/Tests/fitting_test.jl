using JLD2, Optim, CairoMakie, StatsBase, FHist

folder = "525"
file = "1"

parsed_data = JLD2.load_object("../results/"*folder*"/"*file*"/data/parsed_data.jld2")
delays   = parsed_data.delays
arrivals = parsed_data.arrivals

photons_per_bin = 1000
bin_size = 0.001

N = div(length(delays), photons_per_bin)

tau_F = Vector{Float64}(undef, N)
tau_D = Vector{Float64}(undef, N)
A2    = Vector{Float64}(undef, N)
bkg   = Vector{Float64}(undef, N)

Rtime = Vector{Float64}(undef, N)

function MLP(params, time, decay)
    tau_F, tau_D, A2, bkg = params
    if tau_F < 0 || tau_D < 0 || A2 < 0 || bkg < 0 
        return Inf
    end

    L = @. 1/tau_F * exp(-time/tau_F) + A2/tau_D * exp(-time/tau_D) + bkg
    L ./= sum(L)

    return -sum(decay .* log.(L))
end

for i in 1:N
    photons_idx = ((i-1)*photons_per_bin+1) : i*photons_per_bin
    data = delays[photons_idx] * 1e9
    Rtime[i] = arrivals[photons_idx[end]]

    hist_obj = fit(Histogram, data, 0.5 : (ceil(maximum(data))-0.5))
    edges = collect.(hist_obj.edges)[1]

    time = edges[1:end-1] .+ (edges[2] - edges[1])/2
    decay = hist_obj.weights ./ sum(hist_obj.weights)

    x0 = [10, 1e-3, 1, 100]
    results = Optim.optimize(
                             p -> MLP(p, time, decay),
                             x0
                            )
    tau_F[i], tau_D[i], A2[i], bkg[i] = Optim.minimizer(results)
end

nbins = 175
max_tau_F =ceil(maximum(tau_F))
step = max_tau_F / nbins
println(step)

    fig2d = Figure()
    h2d = Hist2D((tau_F, tau_F))
    _, _heatmap = plot(fig2d[1,2], h2d)
    contour!(h2d; levels=[100, 200, 300], labels=true, colormap=:inferno)
    statbox!(fig2d, h2d; position=(1,1))
    Colorbar(fig2d[1,3], _heatmap)
    fig2d

save("hist_tau.png", fig2d)
