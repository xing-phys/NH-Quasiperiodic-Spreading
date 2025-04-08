println("Hello World!")

using CSV, DataFrames


include("../module.jl")
using .basics

F14 = 377
F15 = 610
F16 = 987
F17 = 1597
F18 = 2584
F19 = 4181


L = F16
x0s = ceil(Int, L/2) .+ collect(-50:10:50)
nx = length(x0s)
Nave = 200
n_time = 1000
delta_t = 1

println("I'm happily starting my job!")

V = 0.5 * exp(im * π/3)

@time data = dynamic_exponent_chen(V, Nave, n_time, delta_t, L, x0s)

using Plots

display(scatter(data[:, 1], data[:, 2]))

file_name = "chen_evo.csv"
file_path = "/Users/xingzy/Projects/nh_quasiperiodic_disorder/FigS3/"
CSV.write(file_path * file_name, data)
