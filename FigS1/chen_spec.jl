include("../module.jl")
using .basics

using LinearAlgebra
using Plots

F14 = 377
F15 = 610
F16 = 987
F17 = 1597
F18 = 2584
F19 = 4181

L = F16

V = 0.5 * exp(im * π/3)
ham = basics.chen_obc(V, 10, L, 0.5, 1.5)

e_wf = [eigen(ham[n]) for n in 1:10] 
    
es = vcat([e_wf[n].values for n in 1:10]...)
wfs = hcat([e_wf[n].vectors for n in 1:10]...)

fds = basics.fractal_dim.(eachcol(wfs))

scatter(real.(es), imag.(es), zcolor=fds, markerstrokewidth=0, markersize=2, markeralpha=0.5, clims=(0,1))

using DataFrames, CSV
spec = DataFrame(re=real.(es), im=imag.(es), fd=fds)

file_name = "chen_spec.csv"
file_path = "/Users/xingzy/Projects/nh_quasiperiodic_disorder/FigS3/"
CSV.write(file_path * file_name, spec)
