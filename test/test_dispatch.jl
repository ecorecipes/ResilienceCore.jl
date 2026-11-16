@testset "the right metrics for the right objects" begin
    J = [-1.0 5.0 0.0; 0.0 -2.0 0.0; 0.0 0.0 -0.5]
    A = [0.0 3.0 1.0; 0.5 0.0 0.0; 0.0 0.3 0.0]
    P = [0.9 0.1; 0.2 0.8]
    eq, proj, markov = ContinuousJacobian(J, [1.0, 0.5, 2.0]), ProjectionMatrix(A), MarkovKernel(P)

    # growth belongs to density-independent models; stability and press responses to equilibria
    @test_throws MethodError growth_rate(eq)
    @test_throws MethodError isstable(proj)
    @test_throws MethodError stochastic_invariability(proj)
    @test_throws MethodError press_sensitivity(proj)
    @test_throws MethodError median_reactivity(proj)
    @test_throws MethodError reactivity(markov)
    @test_throws MethodError growth_rate(markov)

    @test growth_rate in applicable_metrics(proj) && !(growth_rate in applicable_metrics(eq))
    @test isstable in applicable_metrics(eq) && press_sensitivity in applicable_metrics(eq)
    @test applicable_metrics(markov) == [resilience, stationary_distribution, mixing_time]

    # the error message names what is defined
    err = try
        growth_rate(eq)
    catch e
        e
    end
    msg = sprint(showerror, err)
    @test occursin("Metrics defined for it", msg) && occursin("resilience", msg)

    # constructors check what the matrix can mean
    @test_throws ArgumentError ProjectionMatrix([1.0 -0.1; 0.2 0.5])
    @test_throws ArgumentError Generator([-1.0 -0.1; 0.2 -0.5])
    @test Generator([-1.0 0.1; 0.2 -0.5]) isa ProjectionOperator{ContinuousTime}
    @test_throws ArgumentError MarkovKernel([0.5 0.4; 0.2 0.8])
    @test_throws ArgumentError ProjectionMatrix(A; weights = [1.0, 0.0, 1.0])
    @test_throws ArgumentError ContinuousJacobian(J, [1.0, 2.0])
    @test_throws ArgumentError median_reactivity(ContinuousJacobian(J))   # no equilibrium, no default pulses
    @test timesemantics(eq) == ContinuousTime() && timesemantics(proj) == DiscreteTime()
end
