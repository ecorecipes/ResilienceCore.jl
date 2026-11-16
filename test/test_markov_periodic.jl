@testset "Markov kernels" begin
    a, b = 0.1, 0.3
    m = MarkovKernel([1-a a; b 1-b])
    @test resilience(m) ≈ -log(abs(1 - a - b))
    @test stationary_distribution(m) ≈ [b, a] ./ (a + b)
    t = mixing_time(m)
    @test t isa Integer && t ≥ 1
    @test mixing_time(m; ε = 1e-12, tmax = 5) === nothing
end

@testset "periodic operators" begin
    A1 = [0.0 3.0 1.0; 0.5 0.0 0.0; 0.0 0.3 0.0]
    A2 = [0.0 1.0 2.0; 0.7 0.0 0.0; 0.0 0.5 0.2]
    p = PeriodicOperator([ProjectionMatrix(A1), ProjectionMatrix(A2)])
    @test period(p) == 2
    @test operator(monodromy(p)) ≈ A2 * A1
    @test growth_rate(p) ≈ sqrt(growth_rate(ProjectionMatrix(A2 * A1)))
    @test resilience(p) ≈ resilience(ProjectionMatrix(A2 * A1)) / 2
    # a constant sequence is the operator itself
    same = PeriodicOperator([ProjectionMatrix(A1), ProjectionMatrix(A1)])
    @test growth_rate(same) ≈ growth_rate(ProjectionMatrix(A1))
    @test resilience(same) ≈ resilience(ProjectionMatrix(A1))
    @test_throws MethodError isstable(p)
    @test_throws MethodError reactivity(p)

    # seasonal fixed point: stable over the period although one season amplifies
    J1 = [1.2 0.0; 0.0 0.5]
    J2 = [0.5 0.0; 0.0 1.1]
    q = PeriodicOperator([DiscreteJacobian(J1), DiscreteJacobian(J2)])
    @test isstable(q) && !isstable(DiscreteJacobian(J1))
    @test resilience(q) ≈ -log(0.6) / 2

    # parts must be of one kind
    @test_throws ArgumentError PeriodicOperator([ProjectionMatrix(A1), DiscreteJacobian(zeros(3, 3))])
end
