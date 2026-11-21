using SciMLBase: SciMLBase, ODEProblem, ODEFunction, SteadyStateProblem, DiscreteProblem
using ForwardDiff: ForwardDiff
using StructuredPopulationCore: ContinuousDomain, DiscreteDomain

@testset "linearise: SciML problems" begin
    # Lotka–Volterra competition: du/dt = r u (1 − A u); equilibrium u* = A⁻¹ 1, J = −diag(r u*) A
    r = [1.0, 0.7, 1.3]
    A = [1.0 0.3 0.2; 0.4 1.0 0.3; 0.1 0.5 1.0]
    ustar = A \ ones(3)
    Jtrue = -Diagonal(r .* ustar) * A
    lv(u, p, t) = p.r .* u .* (1 .- p.A * u)
    lv!(du, u, p, t) = (du .= p.r .* u .* (1 .- p.A * u); nothing)
    p = (; r, A)

    oop = ODEProblem(lv, ustar, (0.0, 10.0), p)
    iip = ODEProblem(lv!, ustar, (0.0, 10.0), p)
    J = linearise(oop)
    @test J isa ContinuousJacobian
    @test J.J ≈ Jtrue
    @test J.equilibrium ≈ ustar
    @test linearise(iip, ustar).J ≈ Jtrue
    @test resilience(J) ≈ -maximum(real, eigvals(Jtrue))
    @test linearise(oop; equilibrium = Float64[]).equilibrium == Float64[]
    @test linearise(SteadyStateProblem(oop)).J ≈ Jtrue

    # an analytic jac is used when present
    calls = Ref(0)
    jac(u, p, t) = (calls[] += 1; -Diagonal(p.r .* u) * p.A + Diagonal(p.r .* (1 .- p.A * u)))
    withjac = ODEProblem(ODEFunction(lv; jac), ustar, (0.0, 10.0), p)
    @test linearise(withjac).J ≈ Jtrue
    @test calls[] == 1

    # not an equilibrium: an informative error, unless the check is switched off
    @test_throws ArgumentError linearise(oop, ustar .+ 0.1)
    @test linearise(oop, ustar .+ 0.1; check = false) isa ContinuousJacobian
    # without ForwardDiff (simulated by an argument it does not handle) the fallback explains what to load
    @test_throws ArgumentError ResilienceCore.autodiff_jacobian(identity, Any[1.0], false)

    # logistic map x′ = a x (1 − x): fixed point 1 − 1/a, multiplier 2 − a
    a = 2.8
    logistic(u, p, t) = p .* u .* (1 .- u)
    logistic!(du, u, p, t) = (du .= p .* u .* (1 .- u); nothing)
    xstar = [1 - 1 / a]
    D = linearise(DiscreteProblem(logistic, xstar, (0.0, 10.0), a))
    @test D isa DiscreteJacobian
    @test D.J ≈ fill(2 - a, 1, 1)
    @test resilience(D) ≈ -log(abs(2 - a))
    @test isstable(D)
    @test linearise(DiscreteProblem(logistic!, xstar, (0.0, 10.0), a)).J ≈ D.J
    @test_throws ArgumentError linearise(DiscreteProblem(logistic, [0.5], (0.0, 10.0), a))

    # a linear flow and its time-Δt map give consistent metrics
    Q = [-0.5 0.4 0.0; 0.3 -0.9 0.2; 0.0 0.5 -0.4]
    Δt = 0.05
    flow = linearise(ODEProblem((u, p, t) -> Q * u, zeros(3), (0.0, 1.0)))
    step = linearise(DiscreteProblem((u, p, t) -> exp(Q * Δt) * u, zeros(3), (0.0, 1.0)))
    @test flow.J ≈ Q
    @test resilience(step) / Δt ≈ resilience(flow) rtol = 1e-8
end

# linearise for matrix and integral projection models is tested in those packages
@testset "mesh_weights" begin
    @test mesh_weights(ContinuousDomain(0.0, 50.0, 100)) ≈ fill(0.5, 100)
    @test mesh_weights(DiscreteDomain([:a, :b])) == ones(2)
end
