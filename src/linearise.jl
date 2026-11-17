# Getting operators out of models. The methods live in package extensions; their documentation
# lives here, on the `linearise` binding, so that it is part of the manual whether or not the
# extensions are loaded.

"""
    linearise(prob::SciMLBase.ODEProblem, u = prob.u0; p = prob.p, t = first(prob.tspan),
              tol = 1e-8, check = true, equilibrium = u)              -> ContinuousJacobian
    linearise(prob::SciMLBase.SteadyStateProblem, u = prob.u0; ...)  -> ContinuousJacobian
    linearise(prob::SciMLBase.DiscreteProblem, u = prob.u0; ...)     -> DiscreteJacobian

The Jacobian of a SciML model at an equilibrium (`f(u) = 0`) or a fixed point of a map
(`f(u) = u`), as an [`EquilibriumJacobian`](@ref). This is the universal route from the ecorecipes
packages: each lowers its models to SciML with `to_ode_problem` or `to_discrete_problem`.

- `f` may be in place or out of place.
- If the problem's function carries an analytic `jac`, it is used. Otherwise the Jacobian is
  computed by automatic differentiation, which needs ForwardDiff to be loaded
  (`using ForwardDiff`; OrdinaryDiffEq loads it).
- `u` must be an equilibrium or fixed point: the residual `‖f(u)‖∞` (or `‖f(u) − u‖∞`) must
  be at most `tol · max(1, ‖u‖∞)`, otherwise an `ArgumentError` is thrown. Find `u` first,
  for example with a `SteadyStateProblem`, or pass `check = false`.
- `equilibrium` is stored with the Jacobian for the median metrics, whose default pulses are
  proportional to it. Pass `Float64[]` when the state is not a vector of abundances.

Available when SciMLBase is loaded (extension `ResilienceCoreSciMLBaseExt`).

---

    linearise(m::MatrixProjectionModels.MatrixProjectionModel)   -> ProjectionMatrix
    linearise(prob::MatrixProjectionModels.MPMProblem)           -> ProjectionMatrix
    linearise(prob::IntegralProjectionModels.IPMProblem)         -> ProjectionMatrix

The density-independent projection operator of a matrix or integral projection model, as a
[`ProjectionOperator`](@ref). For an IPM the kernel is discretised with `materialize` and the
operator carries the mesh weights of its domain (see [`mesh_weights`](@ref)), so that its
metrics approximate those of the integral operator and converge as the mesh is refined.

Only density-independent, deterministic problems have a single projection operator. A
density-dependent model is linearised at a fixed point instead, through its SciML
`DiscreteProblem`; a stochastic one has a stochastic growth rate but no single spectrum.

Available when MatrixProjectionModels or IntegralProjectionModels is loaded (extensions
`ResilienceCoreMatrixProjectionModelsExt` and `ResilienceCoreIntegralProjectionModelsExt`).
"""
linearise

"""
    mesh_weights(domain::ContinuousDomain) -> fill(step_size(domain), n_states(domain))
    mesh_weights(domain::DiscreteDomain)   -> ones(n_states(domain))

The weights of the states of a discretised domain from StructuredPopulationCore: the bin width
`h` for each mesh point of a continuous domain (midpoint rule) and one for each level of a
discrete domain. Passed as `weights` to a [`ProjectionOperator`](@ref), they make norms and
covariances approximate integrals.
"""
mesh_weights(d::ContinuousDomain) = fill(Float64(step_size(d)), n_states(d))
mesh_weights(d::DiscreteDomain) = ones(n_states(d))

# Jacobian of `f` at `u` by automatic differentiation. `inplace = true` means `f(du, u)`;
# otherwise `f(u)` returns the value. Methods with concrete argument types are added by
# ResilienceCoreForwardDiffExt; this fallback explains what to load.
function autodiff_jacobian(f, u, inplace)
    throw(ArgumentError("computing this Jacobian needs automatic differentiation: load ForwardDiff " *
                        "(`using ForwardDiff`; OrdinaryDiffEq loads it), or give the problem's function " *
                        "an analytic `jac`"))
end

"""
    affine_equilibrium(G, b) -> Vector{Float64}

The equilibrium `u* = −G⁻¹ b` of the affine linear model `du/dt = G u + b`, for example a
stage-structured population with constant immigration or recruitment `b`. Its Jacobian is `G`
itself, so `ContinuousJacobian(G, affine_equilibrium(G, b))` describes the return to it.

Throws an `ArgumentError` if `G` is not stable (spectral abscissa `≥ 0`, so no attracting
equilibrium) or is singular. For a generator with non-negative off-diagonal entries that is
stable, `−G⁻¹ ≥ 0`, so a non-negative `b` gives a non-negative `u*`.
"""
function affine_equilibrium(G::AbstractMatrix, b::AbstractVector)
    n = checksquare(G, "G")
    length(b) == n || throw(ArgumentError("b has $(length(b)) entries for a $n × $n matrix"))
    s = spectral_abscissa(G)
    s < 0 || throw(ArgumentError("G is not stable (spectral abscissa $s ≥ 0): du/dt = G u + b has no attracting " *
                                 "equilibrium"))
    F = lu(Matrix{Float64}(G); check = false)
    issuccess(F) || throw(ArgumentError("G is singular: du/dt = G u + b has no unique equilibrium"))
    return -(F \ Vector{Float64}(b))
end

"""
    periodic_generator(G, T; steps = 100, weights = nothing) -> PeriodicOperator

The period map of a density-independent model `du/dt = G(t) u` whose generator `G` (a function of
time returning a matrix with non-negative off-diagonal entries) has period `T`, as a
[`PeriodicOperator`](@ref) of `steps` projection matrices `exp(G(tᵢ) Δt)`, with `Δt = T / steps` and
`tᵢ = (i − ½) Δt` the midpoints of the sub-intervals of `[0, T]`. No ODE solver is needed. The
approximation is exact when the `G(t)` commute (for example `G(t) = f(t) Q`) and the midpoint
rule integrates `f` exactly, and converges as `steps` grows otherwise.

Each part is one sub-step `Δt`, so the operator's per-step metrics are per `Δt`: the growth factor
over a period is `growth_rate(op)^steps`, and the growth rate per unit time is
`log(growth_rate(op)) / Δt`. `weights` are the mesh weights of a discretised continuous state
(see [`mesh_weights`](@ref)); ones by default.
"""
function periodic_generator(G, T::Real; steps::Integer = 100, weights = nothing)
    T > 0 || throw(ArgumentError("the period must be positive, got $T"))
    steps ≥ 1 || throw(ArgumentError("steps must be at least 1, got $steps"))
    Δt = T / steps
    parts = map(1:steps) do i
        t = (i - 0.5) * Δt
        Gt = G(t)
        Q = Generator(Gt; weights = something(weights, ones(size(Gt, 1))))   # checks off-diagonal signs
        E = exp(Δt * operator(Q))
        # exp of a matrix with non-negative off-diagonal entries is non-negative; clear rounding
        tiny = 64 * eps() * max(1.0, maximum(abs, E))
        all(≥(-tiny), E) || throw(ArgumentError("exp(G(t) Δt) has negative entries at t = $t"))
        ProjectionMatrix(max.(E, 0.0); weights = Q.weights)
    end
    return PeriodicOperator(parts)
end
