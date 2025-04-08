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


L = F18
x0s = ceil(Int, L/2) .+ collect(-50:10:50)
nx = length(x0s)
Nave = 50
n_time = 1200
delta_t = .1

println("I'm happily starting my job!")

gamma = π/3
λs = [0.2, 0.5, 1, 1.5, 2.5]

for j in eachindex(λs)
    lambda = λs[j]
    @time data = dynamic_exponent(lambda * exp(im * gamma), Nave, n_time, delta_t, L, x0s)
    file_name = "lambda=$(lambda).csv"
    file_path = "/home/xing/Documents/Julia/nh_quasiperiodic/evolution_data/"
    data = CSV.write(file_path * file_name, data)
    println("λ=$(lambda) done!")
end