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
