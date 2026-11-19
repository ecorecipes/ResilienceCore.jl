# A toy community of independent logistic species, dxᵢ/dt = rᵢxᵢ(1 − xᵢ/Kᵢ) − mᵢdᵢxᵢ, whose
# equilibria and extinctions are known exactly: species i persists while mᵢdᵢ < rᵢ.
struct Logistic
    r::Vector{Float64}
    K::Vector{Float64}
    d::Vector{Float64}
    present::BitVector
    x::Vector{Float64}
end
function Logistic(r, K, d, m, present)
    x = [present[i] ? max(K[i] * (1 - m[i] * d[i] / r[i]), 0.0) : 0.0 for i in eachindex(r)]
    return Logistic(r, K, d, present .& (x .> 0), x)
end
ResilienceCore.state(c::Logistic) = c.x
ResilienceCore.present_species(c::Logistic) = findall(c.present)
ResilienceCore.with_mortality(c::Logistic, f; newton = true, steady = nothing) = Logistic(c.r, c.K, c.d, f, c.present)
function ResilienceCore.without_species(c::Logistic, i; steady = nothing)
    p = copy(c.present)
    p[i] = false
    return Logistic(c.r, c.K, c.d, ones(length(c.r)), p)
end
ResilienceCore.jacobian(c::Logistic) = Diagonal(-(c.r .* c.x ./ c.K)[c.present])

@testset "simulation experiments on a community interface" begin
    r, K, d = [1.0, 1.0, 1.0], [2.0, 1.0, 0.5], [0.099, 0.199, 0.249]
    c = Logistic(r, K, d, ones(3), trues(3))
    @test richness(c) == 3 && abundances(c) == c.x

    # global tolerance: the first extinction is at the smallest multiplier ≥ min rᵢ/dᵢ on the grid
    opts = PerturbationOptions()
    kexp = findfirst(k -> (opts.global_start + opts.step * k) * maximum(d ./ r) ≥ 1, 0:2000) - 1
    @test global_mortality_tolerance(c) ≈ opts.global_start + opts.step * kexp
    # local tolerances: species i dies at the smallest 1 + step k ≥ rᵢ/dᵢ
    loc = local_mortality_tolerances(c)
    @test loc ≈ [opts.step * findfirst(k -> (1 + opts.step * k) * d[i] ≥ r[i], 1:2000) for i in 1:3]

    # a 10% press on everyone: each xᵢ falls by Kᵢ · 0.1 dᵢ/rᵢ
    rel, l1, lost = press_response(c, fill(1.1, 3))
    drop = K .* 0.1 .* d ./ r
    @test rel ≈ sum(drop) / sum(c.x) && l1 ≈ sum(drop) / sum(c.x) && lost == 0

    # removals: independent species, so each removal loses exactly that species
    rem = removal_responses(c)
    @test all(x -> x.survivors == 2, rem)
    @test [x.relative for x in rem] ≈ c.x ./ sum(c.x)
    @test extinction_tolerance(c) == 2.0          # 3 species: two removals leave at most 1.5

    # linearisation: the continuous-time Jacobian at the equilibrium
    op = linearise(c)
    @test op isa ContinuousJacobian && resilience(op) ≈ minimum(r .* c.x ./ K)
end

@testset "the 27 metrics of Domínguez-García et al. (2019) by name" begin
    r, K, d = [1.0, 1.0, 1.0], [2.0, 1.0, 0.5], [0.099, 0.199, 0.249]
    c = Logistic(r, K, d, ones(3), trues(3))
    B = sum(c.x)
    drop = K .* 0.1 .* d ./ r                    # each species' loss under its own 10% press

    g = global_press_response(c)
    @test values(g) == press_response(c, fill(1.1, 3))
    lp = local_press_responses(c)
    @test [x.relative for x in lp] ≈ drop ./ B && [x.l1 for x in lp] ≈ drop ./ B
    @test all(x -> x.extinctions == 0, lp)

    # the Jacobian part: J = −diag(a), so every quantity is known in closed form
    a = r .* c.x ./ K
    jp = jacobian_profile(linearise(c))
    @test keys(jp) == (:R0, :Is, :Amax, :tmax, :Rinf, :MR0, :MAmax, :Mtmax, :MRinf, :S, :mean_sij)
    @test jp.R0 ≈ -minimum(a) && jp.Rinf ≈ minimum(a) && jp.Is ≈ minimum(a)
    @test jp.Amax == 1 && jp.tmax == 0 && jp.MAmax == 1
    @test jp.S ≈ sum(1 ./ a) && jp.mean_sij ≈ sum(1 ./ a) / 9
    @test jp.MR0 ≈ median_reactivity(linearise(c))
    # a map: the same metrics through the discrete-time methods
    jd = jacobian_profile(DiscreteJacobian(exp(Matrix(jacobian(c))), c.x))
    @test jd.Rinf ≈ minimum(a) && jd.R0 ≈ -minimum(a) && jd.Amax == 1

    prof = stability_profile(c)
    @test length(prof) == 27
    @test prof[keys(jp)] == jp
    @test prof.TM == global_mortality_tolerance(c) && prof.min_TM == minimum(local_mortality_tolerances(c))
    @test prof.mean_TM ≈ sum(local_mortality_tolerances(c)) / 3
    @test prof.mean_CE == 0 && prof.max_CE == 0 && prof.TE ≈ 2 / 3   # independent species
    @test prof.RM ≈ sum(drop) / B && prof.SM ≈ sum(drop) / B
    @test prof.mean_RM ≈ sum(drop) / 3B && prof.max_RM ≈ maximum(drop) / B
    @test prof.mean_SM ≈ sum(drop) / 3B && prof.max_SM ≈ maximum(drop) / B
    @test prof.mean_RE ≈ 1 / 3 && prof.max_RE ≈ maximum(c.x) / B
    @test prof.mean_SE ≈ 1 / 3 && prof.max_SE ≈ maximum(c.x) / B
end
