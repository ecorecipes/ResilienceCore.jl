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
