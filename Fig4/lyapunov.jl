using LinearAlgebra
using ProgressMeter
using CSV, DataFrames

include("/home/xing/Projects/nh_quasiperiodic_disorder/module.jl")
using .basics

γ = π/3
V = 2.5 * exp(im * γ)
N= 10000
q = 10


re_vals = collect(-3.01:0.01:3.01)
im_vals = collect(-4.51:0.01:4.51)
re_grid = repeat(re_vals', length(im_vals), 1)
im_grid = repeat(im_vals, 1, length(re_vals))
es = reshape(re_grid .+ im_grid .* im, :)

χe = Vector{Float64}(undef, length(es))
@time @showprogress Threads.@threads for i in eachindex(es)
    χe[i] = Lyapunov_Exponent(V, q, N, es[i])
end


CSV.write("le_2.5.csv", DataFrame(re=real.(es), im=imag.(es), χ=χe))

# data = CSV.read("/home/xing/Projects/nh_quasiperiodic_disorder/Fig3/le_0.5.csv", DataFrame)
# χe = reshape(data[:, 3], size(re_grid))

χe = reshape(χe, size(re_grid))
dos = 1/2π * Laplacian(χe, 0.01, 0.01)
idos = sum(dos, dims=2)[:, 1] * 0.01

CSV.write("idos_2.5.csv", DataFrame(ime=im_vals[2:end-1], idos=idos))

# Draft  ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~

# using Plots

# heatmap(Laplacian(χe, 0.1, 0.1))


# Draft  ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~