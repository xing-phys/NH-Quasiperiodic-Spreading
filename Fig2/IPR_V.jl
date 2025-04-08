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

γs = [0, π/6, π/3, π/2]
λs = collect(0.1:0.1:2)

iprs = Array{Float64}(undef, L, length(λs), length(γs))


for ii in eachindex(γs)
    @showprogress for jj in eachindex(λs)
        γ = γs[ii]
        λ = λs[jj]
        V = λ * exp(im * γ)
        ham = hobc(V, 1, L)
        ψs = eigvecs(ham)
        ipr = IPR.(eachcol(ψs))
        iprs[:, jj, ii] = ipr
    end
end

ave_ipr = sum(iprs, dims = 1)[1, :, :]/L

CSV.write("/home/xing/Documents/Julia/nh_quasiperiodic/fig1/ave_ipr.csv", DataFrame(ave_ipr, :auto))



# Draft
# ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~





# ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~ ~~