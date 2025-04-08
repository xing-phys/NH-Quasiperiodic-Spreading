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


L = 1000
x0s = ceil(Int, L/2) .+ collect(-50:10:50)
nx = length(x0s)
Nave = 100
n_time = 1000
delta_t = 1

println("I'm happily starting my job!")

w = 3

@time data = basics.dynamic_exponent_rand(w, Nave, n_time, delta_t, L, x0s)

using Plots

display(scatter(data[:, 1], data[:, 2]))

file_name = "rand_evo.csv"
file_path = "/Users/xingzy/Projects/nh_quasiperiodic_disorder/FigS3/"
CSV.write(file_path * file_name, data)