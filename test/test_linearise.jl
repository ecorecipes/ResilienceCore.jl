using SciMLBase: SciMLBase, ODEProblem, ODEFunction, SteadyStateProblem, DiscreteProblem
using ForwardDiff: ForwardDiff
using MatrixProjectionModels: MatrixProjectionModels, MatrixProjectionModel, MPMProblem
using IntegralProjectionModels: IntegralProjectionModels, IPMProblem, PKernel, FKernel, MegaKernel, MatrixKernel,
                                LinearSurvival, NormalGrowth, FecundityRate, GeneralIPM, materialize
using StructuredPopulationCore: ContinuousDomain, DiscreteDomain, DensityDependent, StochasticKernelResampled,
                                step_size, meshpoints

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

@testset "linearise: matrix projection models" begin
    A = [0.0 3.0 1.0; 0.5 0.0 0.0; 0.0 0.3 0.0]
    m = MatrixProjectionModel(A)
    @test linearise(m) isa ProjectionMatrix
    @test growth_rate(linearise(m)) ≈ StructuredPopulationCore.lambda(A)
    prob = MPMProblem(m, ones(3), (0, 10))
    @test growth_rate(linearise(prob)) ≈ StructuredPopulationCore.lambda(A)
    @test resilience(linearise(prob)) ≈ log(StructuredPopulationCore.damping_ratio(A))
    @test growth_rate(linearise(MPMProblem(A, ones(3), (0, 10)))) ≈ StructuredPopulationCore.lambda(A)
    @test_throws ArgumentError linearise(MPMProblem(DensityDependent(), (n, t, p) -> A, ones(3), (0, 10)))
    @test_throws ArgumentError linearise(MPMProblem([m, MatrixProjectionModel(0.9A)], ones(3), (0, 10)))
end

@testset "linearise: integral projection models" begin
    s = LinearSurvival(2.2, 0.25)
    g = NormalGrowth(0.2, 1.02, 0.7)
    f = FecundityRate(0.003, 0.015, 2.0, 0.3, 1.0)
    function ipm(n)
        domain = ContinuousDomain(0.0, 50.0, n)
        return IPMProblem(PKernel(s, g, domain) + FKernel(f, domain), domain, ones(n), (0, 10))
    end

    op = linearise(ipm(100))
    @test op isa ProjectionMatrix
    @test op.weights ≈ fill(0.5, 100)
    @test mesh_weights(ContinuousDomain(0.0, 50.0, 100)) ≈ fill(0.5, 100)
    @test mesh_weights(DiscreteDomain([:a, :b])) == ones(2)
    @test growth_rate(op) ≈ StructuredPopulationCore.lambda(materialize(ipm(100).kernel))

    # the metrics approximate those of the integral operator, so they converge with the mesh
    fine = linearise(ipm(400))
    @test growth_rate(op) ≈ growth_rate(fine) rtol = 1e-3
    @test reactivity(op) ≈ reactivity(fine) rtol = 1e-2
    @test max_amplification(op)[1] ≈ max_amplification(fine)[1] rtol = 1e-2

    # a general IPM with a discrete seedbank: weights h for the mesh, 1 for the seedbank
    size_domain = ContinuousDomain(0.0, 10.0, 20)
    states = (size = size_domain, seedbank = DiscreteDomain([:seedbank]))
    h = step_size(size_domain)
    z = meshpoints(size_domain)
    mk = MegaKernel(; states,
                    size_to_size = PKernel(LinearSurvival(1.0, 0.1), NormalGrowth(0.5, 0.9, 0.5), size_domain),
                    seedbank_to_size = MatrixKernel(reshape(0.3 .* exp.(-(z .- 2) .^ 2), 20, 1)),
                    size_to_seedbank = MatrixKernel(reshape(0.1 .* exp.(0.01 .* z) .* h, 1, 20)),
                    seedbank_to_seedbank = MatrixKernel(fill(0.5, 1, 1)))
    general = linearise(IPMProblem(GeneralIPM(), mk, states, ones(21), (0, 10)))
    @test general.weights ≈ [fill(h, 20); 1.0]
    @test growth_rate(general) ≈ StructuredPopulationCore.lambda(materialize(mk))
end
