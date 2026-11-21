# ResilienceCore.jl

Resilience and stability metrics that dispatch on what a matrix **means**.

The same square matrix can be the Jacobian of a flow at an equilibrium, the Jacobian of a map at
a fixed point, a density-independent projection matrix (matrix or integral projection model), a
generator, the period map of a seasonal model, or a Markov kernel. Stability and resilience mean
different things for each. A projection matrix has a growth rate and a damping ratio but no
equilibrium, and a Markov kernel's dominant eigenvalue is always 1. A metric computed on the
wrong kind of object returns a number without meaning.

ResilienceCore puts that meaning in the type. Each metric has methods only for the objects on
which it is defined, so a wrong pairing is a `MethodError` whose message lists the metrics that
*are* defined.

```julia
using ResilienceCore

J = [-1.0 5.0; 0.0 -2.0]
eq = ContinuousJacobian(J, [1.0, 0.5])      # a flow at an equilibrium
resilience(eq), reactivity(eq), max_amplification(eq), median_reactivity(eq)

A = [0.0 3.0 1.0; 0.5 0.0 0.0; 0.0 0.3 0.0]
proj = ProjectionMatrix(A)                   # density-independent, discrete time
growth_rate(proj), resilience(proj)          # λ and the log damping ratio
reactivity(proj)                             # transient growth relative to λ (Stott et al. 2011)

growth_rate(eq)                              # MethodError: ... Metrics defined for it: resilience, isstable, …
applicable_metrics(proj)
```

| object | type | e.g. from |
|---|---|---|
| flow at an equilibrium | `ContinuousJacobian(J, x̄)` | ODE food webs, Lotka–Volterra, epidemic models at an endemic equilibrium |
| map at a fixed point | `DiscreteJacobian(J, x̄)` | density-dependent matrix and integral projection models, matrix community models |
| density-independent projection | `ProjectionMatrix(A; weights)` | MatrixProjectionModels.jl, IntegralProjectionModels.jl (with mesh `weights`) |
| density-independent generator | `Generator(Q; weights)` | FiniteStatePopulationDynamics.jl |
| seasonal model | `PeriodicOperator([A₁, A₂, …])` | periodic matrices, seasonal fixed points |
| Markov kernel | `MarkovKernel(P)` | stochastic state transitions, Bayesian-network kernels |

The time semantics, `ContinuousTime` and `DiscreteTime`, are those of
[StructuredPopulationCore](https://github.com/ecorecipes/ProjectionModels.jl), the shared core
of the ecorecipes projection-model packages. The tests check the identities that tie the
definitions together:
- the metrics of `exp(J Δt)` per `Δt` equal those of `J`;
- transients of a projection matrix are unchanged when it is rescaled;
- growth rate and damping ratio agree with StructuredPopulationCore;
- integral-projection metrics do not depend on the mesh once it is weighted.

## From models to operators

`linearise` builds the operator from a model. The SciML methods come from a package extension;
modelling packages that depend on ResilienceCore add methods for their own problems:

```julia
using ResilienceCore, OrdinaryDiffEq        # loads SciMLBase and ForwardDiff
J = linearise(ODEProblem(f, u_star, tspan, p))           # ContinuousJacobian at u_star
D = linearise(DiscreteProblem(g, x_star, tspan, p))      # DiscreteJacobian at a fixed point

using MatrixProjectionModels, IntegralProjectionModels  # each defines its own methods
linearise(MPMProblem(mpm, n0, (0, 50)))                  # ProjectionMatrix
linearise(IPMProblem(kernel, domain, n0, (0, 50)))       # ProjectionMatrix with mesh weights
```

The SciML route reaches every ecorecipes model with a steady state, because the packages lower
their models with `to_ode_problem` and `to_discrete_problem`. `linearise` checks that the point
is an equilibrium or fixed point, and uses an analytic `jac` when the problem has one.

## Gallery

The [gallery](https://ecorecipes.github.io/ResilienceCore.jl/dev/gallery/) links to the
resilience tutorials of the ecorecipes packages that use ResilienceCore:
- SciML models, through [DynamicalResilience.jl](https://github.com/ecorecipes/DynamicalResilience.jl);
- community matrices, food webs and structured populations;
- Markov chains, sign patterns and agent-based and cellular-automaton models.

## Installation

Neither this package nor StructuredPopulationCore is registered yet:

```julia
using Pkg
Pkg.add(url = "https://github.com/ecorecipes/ProjectionModels.jl")   # StructuredPopulationCore
Pkg.add(url = "https://github.com/ecorecipes/ResilienceCore.jl")
```

## Status

The operator types and their metrics:
- growth and return rates;
- reactivity and amplification in weighted norms;
- the median metrics of Arnoldi et al.;
- press sensitivity;
- stochastic invariability;
- growth sensitivity and elasticity;
- Floquet return rates;
- Markov mixing.

`linearise` for SciML ODE, steady-state and discrete problems, and for matrix and integral
projection models.

Planned, in order:
- extensions for MatrixCommunityModels, FiniteStatePopulationDynamics and
  ContinuousStatePopulationDynamics;
- stochastic environments;
- simulation-based resilience (basin stability, early-warning indicators);
- composition results (Kronecker spectra of stratified models, lumping bounds).
