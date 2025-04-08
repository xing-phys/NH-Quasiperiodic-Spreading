module basics
    using LinearAlgebra, DifferentialEquations, ProgressMeter, DataFrames
    # using ExponentialUtilities
    export hobc, chen_obc, rand_obc, hpbc, hobc_benchmark, IPR, fractal_dim, ComplexRatio, RealRatio, RealStat, DOS, DOS_distribution, evolution_sol, evolution_mat, σ₂, σ, dynamic_exponent, dynamic_exponent_chen, dynamic_exponent_rand, Lyapunov_Exponent, Lyapunov_Exponent_chen, Lyapunov_Exponent_rand, Laplacian, Lyapunov_Evolution

    function hobc(V, Nave, L)
        t = 1 + 0im
        α = (sqrt(5) - 1)/2
        sub = t * ones(L-1)
        sup = t * ones(L-1)
        if Nave == 1
            onsite = 2 * V * cos.(2π * α * collect(1:L))
            ham = Tridiagonal(sub, onsite, sup)
        else
            δ = collect(2π/Nave: 2π/Nave: 2π + 1e-10)
            onsite = [2 * V * cos.(2π * α * collect(1:L) .+ δ[i]) for i in 1:Nave]
            ham = [Tridiagonal(sub, onsite[i], sup) for i in 1:Nave]
        end
        return ham
    end

    function chen_obc(V, Nave, L, b, p)
        t = 1 + 0im
        α = (sqrt(5) - 1)/2
        sub = t * ones(L-1)
        sup = t * ones(L-1)
        if Nave == 1
            onsite = (2 * V * cos.(2π * α * collect(1:L)) .+ p) ./ (1 .- b * cos.(2π * α * collect(1:L)))
            ham = Tridiagonal(sub, onsite, sup)
        else
            δ = collect(2π/Nave: 2π/Nave: 2π + 1e-10)
            onsite = [(2 * V * cos.(2π * α * collect(1:L) .+ δ[i]) .+ p) ./ (1 .- b * cos.(2π * α * collect(1:L) .+ δ[i])) for i in 1:Nave]
            ham = [Tridiagonal(sub, onsite[i], sup) for i in 1:Nave]
        end
        return ham
    end

    function rand_obc(w, L)
        t = 1 + 0im
        sub = t * ones(L-1)
        sup = t * ones(L-1)
        r = w .* sqrt.(rand(L))
        ϕ = 2π .* rand(L)
        onsite = r .* exp.(im .* ϕ)
        ham = Tridiagonal(sub, onsite, sup)
        return ham
    end

    function hpbc(V, Nave, L)
        t = 1 + 0im
        α = (sqrt(5) - 1)/2
        δ = collect(2π/Nave: 2π/Nave: 2π + 1e-10)
        sub = t * ones(L-1)
        sup = t * ones(L-1)
        onsite = [2 * V * cos.(2π * α * collect(1:L) .+ δ[i]) for i in 1:Nave]
        pbc = zeros(ComplexF64, L, L)
        pbc[1, end] = pbc[end, 1] = t
        ham = [Tridiagonal(sub, onsite[i], sup) + pbc for i in 1:Nave]
        return ham
    end

    function hobc_benchmark(V, Nave, L)
        J = 1
        α = (sqrt(5) - 1)/2
        δ = collect(2π/Nave: 2π/Nave: 2π + 1e-10)
        sub = J * ones(L-1)
        sup = J * ones(L-1)
        onsite = [2 * V * cos.(2π * α * collect(1:L) .+ δ[i]) for i in 1:Nave]
        ham = [Tridiagonal(sub, onsite[i], sup) for i in 1:Nave]
        return ham
    end

    function IPR(ψ)
        return sum(x -> abs2(x)^2, ψ)
    end

    function fractal_dim(ψ)
        n = length(ψ)
        return - log(IPR(ψ)) / log(n)
    end

    function ComplexRatio(es)
        N = length(es)
        ratios = Vector{ComplexF64}(undef, N)
        @showprogress Threads.@threads for ii in eachindex(es)
            edif = es .- es[ii]
            sort!(edif, by = abs2)
            ratios[ii] = edif[2] / edif[3]
        end
        re_vals = collect(-1.0:0.02:1.0)
        im_vals = collect(-1.0:0.02:1.0)
        re_grid = repeat(re_vals', length(im_vals), 1)
        im_grid = repeat(im_vals, 1, length(re_vals))
        complex_grid = reshape(re_grid .+ im_grid .* im, :)
        density = Vector{Float64}(undef, length(complex_grid))
        @showprogress Threads.@threads for jj in eachindex(complex_grid)
            zz = complex_grid[jj]
            num = count(x -> real(zz)-0.005 < real(x) < real(zz)+0.005 && imag(zz)-0.005 < imag(x) < imag(zz)+0.005, ratios)
            density[jj] = num /N
        end
        density_df = DataFrame(Re=real.(complex_grid), Im=imag.(complex_grid), den=density)
        return density_df
    end

    function RealRatio(es)
        N = length(es)
        ratios = Vector{Float64}(undef, N)
        @showprogress Threads.@threads for ii in eachindex(es)
            edif = abs.(es .- es[ii])
            sort!(edif)
            ratios[ii] = edif[2] / edif[3]
        end
        re_vals = collect(0:0.001:1.0)
        density = Vector{Float64}(undef, length(re_vals))
        @showprogress Threads.@threads for jj in eachindex(re_vals)
            zz = re_vals[jj]
            num = count(x -> zz-0.0005 < x < zz+0.0005, ratios)
            density[jj] = num / N / 0.001
        end
        density_df = DataFrame(Re=re_vals, den=density)
        return density_df
    end

    function RealStat(es, intvl)
        N = length(es)
        sortedes = sort(es)
        ΔE = sortedes[2:end] .- sortedes[1:end-1]
        aveΔE = sum(ΔE) / (N - 1)
        δE = ΔE / aveΔE
        re_vals = collect(0:0.01:intvl)
        density = Vector{Float64}(undef, length(re_vals))
        @showprogress Threads.@threads for jj in eachindex(re_vals)
            zz = re_vals[jj]
            num = count(x -> zz-0.005 < x < zz+0.005, δE)
            density[jj] = num / (N-1) / 0.01
        end
        density_df = DataFrame(Re=re_vals, den=density)
        return density_df
        # return δE
    end

    function DOS(e, de, es)
        # return the dos at energy e with window size de
        return count(x -> e - de/2 < x < e + de/2, es) / length(es)
    end

    function DOS_distribution(es, parts)
        # return the dos and the corresponding energies
        emax = maximum(es)
        emin = minimum(es)
        Δe = (emax - emin) / parts
        start_e = emin + Δe/2
        final_e = emax - Δe/2
        sepes = collect(start_e:Δe:final_e)
        return map(x -> DOS(x, Δe, es), sepes), sepes
    end

    function evolution_sol(ham, ψ, n_time, δt)
        A = -im * ham
        function schrodinger!(dx, x, p, t)
            dx .= A * x
        end
        L = length(ψ)
        ψt = copy(ψ)
        ψts = Matrix{ComplexF64}(undef, L, n_time)
        for j in 1:n_time
            tspan = (0.0, δt)
            sol = solve(ODEProblem(schrodinger!, ψt, tspan))
            tt = sol.t
            nt = length(tt)
            solψ = hcat(sol.u...)
            for j in 1:nt
                solψ[:, j] = normalize(solψ[:, j])
            end
            ψt = normalize(sol.u[end])
            ψts[:, j] = ψt
        end
        return ψts
    end

    function evolution_mat(ham, ψ, n_time, δt)
        L = length(ψ)
        A = -im * ham
        ts = collect(1: n_time)
        ψt = copy(ψ)
        ψts = Matrix{ComplexF64}(undef, L, n_time)
        for t in ts
            ψt = expv(δt, A, ψt)
            normalize!(ψt)
            ψts[:, t] = ψt
        end
        # ψt2 = abs2.(ψts)
        return ψts
    end

    function σ₂(ψ2, x0)
        l = length(ψ2)
        return sum([(i - x0)^2 * ψ2[i] for i in 1:l])
    end

    function σ(ψ2, x0)
        σ2 = σ₂(ψ2, x0)
        return sqrt(σ2)
    end

    function dynamic_exponent(V, Nave, n_time, δt, L, x0s)
        # This function will generate a datafram with log(t) and log(σ)

        # V = λ * exp(im * γ)
        ham = hobc(V, Nave, L)
        # Set the max evolution time
        tmax = n_time * δt
        ts = collect(range(0, tmax, n_time+1))

        nx = length(x0s)

        loop_kernel = [(i, j) for i in 1:Nave, j in 1:nx]
        
        ψt2 = Array{Float64}(undef, L, n_time, Nave, nx)
        @time Threads.@threads for p in eachindex(loop_kernel)
            lk = loop_kernel[p]
            n = lk[1]
            x = lk[2]
            x0 = x0s[x]
            init_state = 0im * zeros(L)
            init_state[x0] = 1
            ψtnx = evolution_sol(ham[n], init_state, n_time, δt)
            ψt2[:, :, n, x] = abs2.(ψtnx)
        end
        ψt2mean = map(ψs -> sum(ψs, dims=3)[:,:,1]/Nave, eachslice(ψt2, dims=4))
        σs = Matrix{Float64}(undef, n_time, nx)

        println("Begin to calclulate mean σs")
        for i in eachindex(x0s)
            σs[:, i] = σ.(eachcol(ψt2mean[i]), x0s[i])
        end
        σmean = sum(σs, dims=2)[:, 1] / nx

        println("Begin to generate data")
        data = DataFrame(X = log.(ts[2:end]), Y = log.(σmean))

        return data
    end

    function dynamic_exponent_chen(V, Nave, n_time, δt, L, x0s)
        # This function will generate a datafram with log(t) and log(σ)

        # V = λ * exp(im * γ)
        ham = chen_obc(V, Nave, L, 0.5, 1.5)
        # Set the max evolution time
        tmax = n_time * δt
        ts = collect(range(0, tmax, n_time+1))

        nx = length(x0s)

        loop_kernel = [(i, j) for i in 1:Nave, j in 1:nx]
        
        ψt2 = Array{Float64}(undef, L, n_time, Nave, nx)
        @time @showprogress Threads.@threads for p in eachindex(loop_kernel)
            lk = loop_kernel[p]
            n = lk[1]
            x = lk[2]
            x0 = x0s[x]
            init_state = 0im * zeros(L)
            init_state[x0] = 1
            ψtnx = evolution_sol(ham[n], init_state, n_time, δt)
            ψt2[:, :, n, x] = abs2.(ψtnx)
        end
        ψt2mean = map(ψs -> sum(ψs, dims=3)[:,:,1]/Nave, eachslice(ψt2, dims=4))
        σs = Matrix{Float64}(undef, n_time, nx)

        println("Begin to calclulate mean σs")
        for i in eachindex(x0s)
            σs[:, i] = σ.(eachcol(ψt2mean[i]), x0s[i])
        end
        σmean = sum(σs, dims=2)[:, 1] / nx

        println("Begin to generate data")
        data = DataFrame(X = log.(ts[2:end]), Y = log.(σmean))

        return data
    end

    function dynamic_exponent_rand(w, Nave, n_time, δt, L, x0s)
        # This function will generate a datafram with log(t) and log(σ)

        # V = λ * exp(im * γ)
        # Set the max evolution time
        tmax = n_time * δt
        ts = collect(range(0, tmax, n_time+1))

        nx = length(x0s)

        loop_kernel = [(i, j) for i in 1:Nave, j in 1:nx]
        
        ψt2 = Array{Float64}(undef, L, n_time, Nave, nx)
        @time @showprogress Threads.@threads for p in eachindex(loop_kernel)
            lk = loop_kernel[p]
            n = lk[1]
            x = lk[2]
            x0 = x0s[x]
            init_state = 0im * zeros(L)
            init_state[x0] = 1
            ham = rand_obc(w, L)
            ψtnx = evolution_sol(ham, init_state, n_time, δt)
            ψt2[:, :, n, x] = abs2.(ψtnx)
        end
        ψt2mean = map(ψs -> sum(ψs, dims=3)[:,:,1]/Nave, eachslice(ψt2, dims=4))
        σs = Matrix{Float64}(undef, n_time, nx)

        println("Begin to calclulate mean σs")
        for i in eachindex(x0s)
            σs[:, i] = σ.(eachcol(ψt2mean[i]), x0s[i])
        end
        σmean = sum(σs, dims=2)[:, 1] / nx

        println("Begin to generate data")
        data = DataFrame(X = log.(ts[2:end]), Y = log.(σmean))

        return data
    end

    function Lyapunov_Exponent(V::ComplexF64, q::Int64, N::Int64, E::ComplexF64)
        # This will return the Lyapunov exponent χ(E)
        L = q * N
        α = (sqrt(5) - 1)/2
        Vns = 2 * V * cos.(2π * α * collect(1:L))
        Q = [1 0; 0 1]
        γ1 = 0
        γ2 = 0
        for i in 1:N
            Tns = [[E - Vns[i*q + 1 - j] -1; 1 0] for j in 1:q]
            Tqs = reduce(*, Tns) * Q
            Tqr = qr(Tqs)
            Q = Tqr.Q
            γs = log.(abs.(diag(Tqr.R)))
            γ1 += γs[1]
            γ2 += γs[2]
        end
        γ1 = γ1 / L
        γ2 = γ2 / L
        χ = max(γ1, γ2)
        return χ
    end

    function Lyapunov_Exponent_chen(V::ComplexF64, q::Int64, N::Int64, E::ComplexF64)
        # This will return the Lyapunov exponent χ(E)
        p = 1.5
        b = 0.5
        L = q * N
        α = (sqrt(5) - 1)/2
        Vns = (2 * V * cos.(2π * α * collect(1:L)) .+ p) ./ (1 .- b * cos.(2π * α * collect(1:L)))
        Q = [1 0; 0 1]
        γ1 = 0
        γ2 = 0
        for i in 1:N
            Tns = [[E - Vns[i*q + 1 - j] -1; 1 0] for j in 1:q]
            Tqs = reduce(*, Tns) * Q
            Tqr = qr(Tqs)
            Q = Tqr.Q
            γs = log.(abs.(diag(Tqr.R)))
            γ1 += γs[1]
            γ2 += γs[2]
        end
        γ1 = γ1 / L
        γ2 = γ2 / L
        χ = max(γ1, γ2)
        return χ
    end

    function Lyapunov_Exponent_rand(w::Float64, q::Int64, N::Int64, E::ComplexF64)
        # This will return the Lyapunov exponent χ(E)
        L = q * N
        r = w .* sqrt.(rand(L))
        ϕ = 2π .* rand(L)
        Vns = r .* exp.(im .* ϕ)
        Q = [1 0; 0 1]
        γ1 = 0
        γ2 = 0
        for i in 1:N
            Tns = [[E - Vns[i*q + 1 - j] -1; 1 0] for j in 1:q]
            Tqs = reduce(*, Tns) * Q
            Tqr = qr(Tqs)
            Q = Tqr.Q
            γs = log.(abs.(diag(Tqr.R)))
            γ1 += γs[1]
            γ2 += γs[2]
        end
        γ1 = γ1 / L
        γ2 = γ2 / L
        χ = max(γ1, γ2)
        return χ
    end

    function Laplacian(ϕ, dre, dim)
        N_im, N_re = size(ϕ)
        dos = Matrix{Float64}(undef, N_im-2, N_re-2)
        @time @showprogress Threads.@threads for i in 2:N_re-1
            for j in 2:N_im-1
                ϕ0 = ϕ[j, i]
                ϕ1 = ϕ[j, i - 1]
                ϕ2 = ϕ[j, i + 1]
                ϕ3 = ϕ[j + 1, i]
                ϕ4 = ϕ[j - 1, i]
                dos[j-1, i-1] = (ϕ1 + ϕ2 - 2*ϕ0) / dre^2 + (ϕ3 + ϕ4 - 2*ϕ0) / dim^2
            end
        end
        return dos
    end

    function Lyapunov_Evolution(E, V, N, ϕ, r0)
        α = (sqrt(5) - 1)/2
        rs = Vector{ComplexF64}(undef, N)
        Vs = 2 * V * cos.(2π * α * collect(1:N) .+ ϕ)
        rs[1] = (E - Vs[1]) - 1/r0
        for n in 2:N
            rs[n] = (E - Vs[n]) - 1/rs[n-1]
        end
        return rs
    end
end