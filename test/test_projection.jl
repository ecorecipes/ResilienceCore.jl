@testset "projection matrices and generators" begin
    A = [0.0 3.0 1.0; 0.5 0.0 0.0; 0.0 0.3 0.0]
    proj = ProjectionMatrix(A)

    # agree with StructuredPopulationCore, the shared core of the projection-model packages
    @test growth_rate(proj) ≈ StructuredPopulationCore.lambda(A)
    @test resilience(proj) ≈ log(StructuredPopulationCore.damping_ratio(A))

    # transients are relative to growth: scaling A changes λ, not the transients
    scaled = ProjectionMatrix(2.5A)
    @test growth_rate(scaled) ≈ 2.5growth_rate(proj)
    @test reactivity(scaled) ≈ reactivity(proj)
    @test max_amplification(scaled)[1] ≈ max_amplification(proj)[1]
    @test resilience(scaled) ≈ resilience(proj)

    # the amplification envelope settles at the norm of the dominant-mode projector, not at 0
    Amax, k = max_amplification(proj)
    @test Amax ≥ amplification_envelope(proj, 200) > 0

    # elasticities sum to one
    @test sum(elasticity(proj)) ≈ 1
    @test growth_sensitivity(proj)[2, 1] > 0

    # a generator and its time-Δt projection matrix
    Q = [-0.5 0.4 0.0; 0.3 -0.6 0.2; 0.0 0.5 -0.1]
    gen = Generator(Q)
    Δt = 0.05
    step = ProjectionMatrix(exp(Q * Δt))
    @test log(growth_rate(step)) / Δt ≈ growth_rate(gen) rtol = 1e-10
    @test resilience(step) / Δt ≈ resilience(gen) rtol = 1e-8
    @test reactivity(step) / Δt ≈ reactivity(gen) rtol = 5e-2
end

@testset "integral projection kernels: weighted norms make metrics mesh-independent" begin
    φ(z) = exp(-z^2 / 2) / sqrt(2π)
    K(y, x) = (0.3 + 0.6x) * φ((y - 0.2 - 0.7x) / 0.1) / 0.1 + 2x^2 * φ((y - 0.1) / 0.05) / 0.05
    function discretise(edges)
        mid = (edges[1:(end - 1)] .+ edges[2:end]) ./ 2
        h = diff(edges)
        A = [K(y, x) * hx for y in mid, (x, hx) in zip(mid, h)]   # n(yᵢ) ≈ Σⱼ K(yᵢ, xⱼ) hⱼ n(xⱼ)
        return A, h
    end
    Au, hu = discretise(range(0, 1; length = 401))
    An, hn = discretise(range(0, 1; length = 401) .^ 2)              # mesh concentrated near 0
    uniform = ProjectionMatrix(Au; weights = hu)
    nonuniform = ProjectionMatrix(An; weights = hn)
    @test growth_rate(nonuniform) ≈ growth_rate(uniform) rtol = 1e-2
    @test reactivity(nonuniform) ≈ reactivity(uniform) rtol = 2e-2
    @test max_amplification(nonuniform)[1] ≈ max_amplification(uniform)[1] rtol = 2e-2
    # without the mesh weights the population size is miscounted and the metric moves
    @test !isapprox(reactivity(ProjectionMatrix(An)), reactivity(uniform); rtol = 5e-2)
end
