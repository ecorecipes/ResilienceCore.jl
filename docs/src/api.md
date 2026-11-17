# API Reference

## Operators

```@docs
LinearisedDynamics
EquilibriumJacobian
ProjectionOperator
PeriodicOperator
MarkovKernel
timesemantics
operator
period
monodromy
```

## Asymptotic metrics

```@docs
growth_rate
resilience
isstable
```

## Transient metrics

```@docs
reactivity
propagator
amplification_envelope
max_amplification
```

## Median metrics

```@docs
pulse_covariance
median_displacement
median_reactivity
median_max_amplification
median_resilience
```

## Press and stochastic perturbations

```@docs
sensitivity_matrix
press_sensitivity
total_sensitivity
stochastic_invariability
growth_sensitivity
elasticity
```

## Markov kernels

```@docs
stationary_distribution
mixing_time
```

## Community interface and simulation experiments

A model package implements `state`, `present_species`, `with_mortality`, `without_species` and
`jacobian` for its community type; the experiments then run on it unchanged.

```@docs
state
present_species
with_mortality
without_species
jacobian
richness
abundances
PerturbationOptions
global_mortality_tolerance
local_mortality_tolerances
press_response
removal_responses
extinction_tolerance
```

## Operators from models

`linearise` turns a model into the operator its metrics need. A community gives the Jacobian at
its current state. A SciML `ODEProblem`, `SteadyStateProblem` or `DiscreteProblem` gives the
Jacobian at an equilibrium or fixed point. Every ecorecipes package lowers its models to these
with `to_ode_problem` or `to_discrete_problem`, so this route reaches any model with a steady
state. Density-independent matrix and integral projection models give a projection operator,
with mesh weights for an IPM. The model-specific methods are package extensions, loaded with
SciMLBase (plus ForwardDiff for automatic differentiation), MatrixProjectionModels or
IntegralProjectionModels.

```julia
using ResilienceCore, OrdinaryDiffEq          # OrdinaryDiffEq loads SciMLBase and ForwardDiff

lv(u, p, t) = p.r .* u .* (1 .- p.A * u)      # Lotka–Volterra competition
p = (r = [1.0, 0.7], A = [1.0 0.3; 0.4 1.0])
ustar = p.A \ ones(2)
J = linearise(ODEProblem(lv, ustar, (0.0, 1.0), p))   # a ContinuousJacobian
resilience(J), reactivity(J)
```

```@docs
linearise
mesh_weights
```

## Applicability

```@docs
applicable_metrics
```
