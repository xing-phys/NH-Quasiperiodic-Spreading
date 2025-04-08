using LinearAlgebra
using Plots

include("../module.jl")
using .basics

L = 1000
w = 3


es = Vector{ComplexF64}(undef, 10 * L)
wfs = Matrix{ComplexF64}(undef, L, 10 * L)
for ii in 1:10
    ham = basics.rand_obc(w, L)
    e_wf = eigen(ham)
    es[L*(ii-1)+1:L*ii] = e_wf.values
    wfs[:, L*(ii-1)+1:L*ii] = e_wf.vectors
end

fds = basics.fractal_dim.(eachcol(wfs))

scatter(real.(es), imag.(es), zcolor=fds, markerstrokewidth=0, markersize=2, markeralpha=0.5, clims=(0,1))

using DataFrames, CSV
spec = DataFrame(re=real.(es), im=imag.(es), fd=fds)

file_name = "rand_spec.csv"
file_path = "/Users/xingzy/Projects/nh_quasiperiodic_disorder/FigS3/"
CSV.write(file_path * file_name, spec)

