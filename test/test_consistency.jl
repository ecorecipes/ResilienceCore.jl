@testset "a flow and its time-Δt map give the same metrics" begin
    # a stable, non-normal (reactive) Jacobian
    J = [-1.0 5.0 0.0; 0.0 -2.0 0.0; 0.0 0.0 -0.5]
    B = [1.0, 0.5, 2.0]
    C = pulse_covariance(B)
    flow = ContinuousJacobian(J, B)
    @test isstable(flow) && resilience(flow) ≈ 0.5
    @test reactivity(flow) > 0

    # small Δt: rates per unit time agree
    Δt = 1e-3
    map = DiscreteJacobian(exp(J * Δt), B)
    @test resilience(map) / Δt ≈ resilience(flow) rtol = 1e-10        # exact: |e^{λΔt}| = e^{Re λ Δt}
    @test reactivity(map) / Δt ≈ reactivity(flow) rtol = 1e-2
    @test median_reactivity(map, C) / Δt ≈ median_reactivity(flow, C) rtol = 1e-2
    @test stochastic_invariability(map) / Δt ≈ stochastic_invariability(flow) rtol = 1e-2

    # the same time grid: amplification and median metrics agree exactly
    dt = 0.1
    map01 = DiscreteJacobian(exp(J * dt), B)
    Af, tf = max_amplification(flow; dt = dt)
    Am, km = max_amplification(map01)
    @test Af ≈ Am && tf ≈ km * dt
    Mf, Mtf = median_max_amplification(flow, C; dt = dt)
    Mm, Mkm = median_max_amplification(map01, C)
    @test Mf ≈ Mm && Mtf ≈ Mkm * dt
    @test median_resilience(flow, C; dt = dt) ≈ median_resilience(map01, C) / dt rtol = 1e-10

    # press responses: the Euler map x′ = x + Δt f(x) has the same equilibrium response
    euler = DiscreteJacobian(I + Δt * J)
    @test sensitivity_matrix(euler) * Δt ≈ sensitivity_matrix(flow)
    @test press_sensitivity(euler) * Δt ≈ press_sensitivity(flow)

    # J = −a I: not reactive, resilience a, invariability a, sensitivity I/a
    a = 0.7
    iso = ContinuousJacobian(-a * Matrix(1.0I, 4, 4))
    @test resilience(iso) ≈ a && reactivity(iso) ≈ -a
    @test max_amplification(iso) == (1.0, 0.0)
    @test stochastic_invariability(iso) ≈ a
    @test press_sensitivity(iso) ≈ 4 / a

    # dense and iterative invariability agree
    R = [-3.0 1.0 0.5 0.0; 0.2 -2.0 0.0 1.0; 0.0 0.3 -4.0 0.2; 0.1 0.0 0.4 -2.5]
    @test stochastic_invariability(ContinuousJacobian(R); dense_max = 0) ≈
          stochastic_invariability(ContinuousJacobian(R)) rtol = 1e-6
end
