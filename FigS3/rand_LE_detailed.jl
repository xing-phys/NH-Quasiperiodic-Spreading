using LinearAlgebra
using ProgressMeter
using CSV, DataFrames

include("/home/xing/Projects/nh_quasiperiodic_disorder/module.jl")
using .basics

w = 3.
N= 10000
q = 10


re_vals = collect(-3.001:0.001:3.001)
im_vals = collect(1.999:0.001:3.001)
re_grid = repeat(re_vals', length(im_vals), 1)
im_grid = repeat(im_vals, 1, length(re_vals))
es = reshape(re_grid .+ im_grid .* im, :)

χe = Vector{Float64}(undef, length(es))
@time @showprogress Threads.@threads for i in eachindex(es)
    χe[i] = Lyapunov_Exponent_rand(w, q, N, es[i])
end


# CSV.write("le_rand_detailed.csv", DataFrame(re=real.(es), im=imag.(es), χ=χe))

# data = CSV.read("/home/xing/Projects/nh_quasiperiodic_disorder/Fig3/le_0.5.csv", DataFrame)
# χe = reshape(data[:, 3], size(re_grid))

χe = reshape(χe, size(re_grid))
dos = 1/2π * Laplacian(χe, 0.001, 0.001)
idos = sum(dos, dims=2)[:, 1] * 0.001

CSV.write("idos_rand_detailed.csv", DataFrame(ime=im_vals[2:end-1], idos=idos))