include("../module.jl")
using .basics

using LinearAlgebra
using CSV, DataFrames
using ProgressMeter

F14 = 377
F15 = 610
F16 = 987
F17 = 1597
F18 = 2584
F19 = 4181

L = F18

re_V = collect(0:0.05:2.5)
im_V = collect(0:0.05:2.5)
re_grid = repeat(re_V', length(im_V), 1)
im_grid = repeat(im_V, 1, length(re_V))
Vs = reshape(re_grid .+ im_grid .* im, :)

iprs = Matrix{Float64}(undef, L, length(im_V)* length(re_V))

@showprogress for ii in eachindex(Vs)
    V = Vs[ii]
    ham = hobc(V, 1, L)
    ψs = eigvecs(ham)
    ipr = IPR.(eachcol(ψs))
    iprs[:, ii] = ipr
end

ave_ipr = sum(iprs, dims=1)[1, :]/L

CSV.write("phase_ipr.csv", DataFrame("Re" => real.(Vs), "Im" => imag.(Vs), "IPR" => ave_ipr))

println("Congratulations!")