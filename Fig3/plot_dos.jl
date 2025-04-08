using Plots
using CSV, DataFrames

include("/home/xing/Projects/nh_quasiperiodic_disorder/module.jl")
using .basics

data = CSV.read("/home/xing/Projects/nh_quasiperiodic_disorder/Fig3/le_0.5.csv", DataFrame)

χe = Matrix(data)
res = unique(χe[:, 1])
ims = unique(χe[:, 2])
χs = reshape(χe[:, 3], length(ims), length(res))

dos = 1/2π * Laplacian(χs[1:2:end, 1:2:end], 0.02, 0.02)


plot_χe = heatmap(res,ims, χs, aspect_ratio=:equal, xlims=(res[1], res[end]), color=cgrad([RGBA(0,1,0,1), RGBA(0,0,1,1)]), title="Lyapunov Exponent χ(E)")




plot_dos = heatmap(res[3:2:end-2], ims[3:2:end-2], dos[1:1:end, 1:1:end], clims=(0, Inf), xlims=(res[2], res[end-1]), color=cgrad([RGB(.9,.9,1), RGB(0,0,0)], [0.0, 0.1]), legend=:false, framestyle=:none, size=(300, 300), background=:transparent)

savefig(plot_dos, "./Fig3/dos_0.5.svg")
