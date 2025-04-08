println("Hello World!")

using ProgressMeter

using LinearAlgebra
using CSV, DataFrames


# ------------------------------------------------------
# Using for linear regression and plot
# using GLM
# using Plots
# ------------------------------------------------------
include("../module.jl")
using .basics

F14 = 377
F15 = 610
F16 = 987
F17 = 1597
F18 = 2584
F19 = 4181


L = F18

println("I'm happily starting my job!")

re_V = collect(0:0.1:2.5)
im_V = collect(0:0.1:2.5)
re_grid = repeat(re_V', length(im_V), 1)
im_grid = repeat(im_V, 1, length(re_V))
Vsall = reshape(re_grid .+ im_grid .* im, :)

Vs = Vsall[1:10]


@showprogress for i in eachindex(Vs)
    V = Vs[i]
    if abs2(V) > 1
        x0s = ceil(Int, L/2) .+ collect(-50:2:50)
        nx = length(x0s)
        Nave = 50
        n_time = 2000
        delta_t = 1.0
    else
        x0s = ceil(Int, L/2) .+ collect(-50:5:50)
        nx = length(x0s)
        Nave = 50
        n_time = 4000
        delta_t = .5
    end
    df = dynamic_exponent(V, Nave, n_time, delta_t, L, x0s)
    file_name = "re$(real(V))_im$(imag(V)).csv"
    CSV.write(file_name, df)
end


# ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~

println("Congratulations!")



# ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~
# Draft
# ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~


# ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~