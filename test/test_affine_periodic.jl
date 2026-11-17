@testset "affine_equilibrium" begin
    G = [-1.0 0.5; 0.2 -0.8]
    b = [1.0, 0.5]
    u = affine_equilibrium(G, b)
    @test G * u + b ≈ zeros(2) atol = 1e-12
    @test all(≥(0), u)
    J = ContinuousJacobian(G, u)
    @test isstable(J)
    @test resilience(J) ≈ -maximum(real, eigvals(G))
    @test_throws ArgumentError affine_equilibrium([0.1 0.0; 0.0 -1.0], b)    # unstable
    @test_throws ArgumentError affine_equilibrium([-1.0 0.0; 0.0 0.0], b)    # singular: abscissa 0
    @test_throws ArgumentError affine_equilibrium(G, [1.0])
end

@testset "periodic_generator" begin
    Q = [-1.0 0.3; 0.5 -0.4]
    T = 2.0
    f(t) = 1 + 0.5 * sin(π * t / T)                     # ∫₀ᵀ f = T + T/π
    exact = exp(Q * (T + T / π))                         # the G(t) = f(t) Q commute
    err(k) = norm(operator(monodromy(periodic_generator(t -> f(t) * Q, T; steps = k))) - exact)
    e = err.((10, 20, 40, 80))
    @test all(e[i + 1] < e[i] for i in 1:3)              # converges
    @test e[end] / e[end - 1] ≈ 0.25 atol = 0.02          # midpoint rule: second order
    @test err(400) < 1e-5

    op = periodic_generator(t -> f(t) * Q, T; steps = 50)
    @test op isa PeriodicOperator
    @test period(op) == 50
    Δt = T / 50
    @test log(growth_rate(op)) / Δt ≈ maximum(real, eigvals(Q)) * (T + T / π) / T rtol = 1e-3

    # non-commuting generators converge too
    R(t) = [-1.0-sin(2π * t / T) 0.3; 0.5+0.4cos(2π * t / T) -0.4]
    M(k) = operator(monodromy(periodic_generator(R, T; steps = k)))
    @test norm(M(400) - M(200)) < 1e-4

    w = [0.5, 0.5]
    @test all(p -> p.weights == w, periodic_generator(t -> Q, T; steps = 3, weights = w).parts)
    @test_throws ArgumentError periodic_generator(t -> [-1.0 -0.1; 0.0 -1.0], T)   # not a generator
    @test_throws ArgumentError periodic_generator(t -> Q, -1.0)
end
