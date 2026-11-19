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
global_press_response
local_press_responses
```

## The 27 metrics by name

See the [Metric catalogue](catalogue.md).

```@docs
jacobian_profile
stability_profile
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

Two helpers serve linear continuous-time models that are not plain generators: an affine
model `du/dt = G u + b` returns to the equilibrium [`affine_equilibrium`](@ref) gives, with
Jacobian `G`, and a periodically forced generator `G(t)` becomes a [`PeriodicOperator`](@ref)
with [`periodic_generator`](@ref), without an ODE solver. Model packages call them from their
own `linearise` methods.

```julia
G = [-1.0 0.5; 0.2 -0.8]
J = ContinuousJacobian(G, affine_equilibrium(G, [1.0, 0.5]))
op = periodic_generator(t -> (1 + 0.5sin(π * t)) * G, 2.0; steps = 100)
```

```@docs
linearise
mesh_weights
affine_equilibrium
periodic_generator
```

## Sign patterns

When only the signs of a Jacobian's entries are known, for example from a causal loop diagram
or a qualitative food-web model, a [`SignPattern`](@ref) supports qualitative stability and
press analysis by sampling magnitudes. Links of unknown sign are sampled with either sign.

```julia
s = SignPattern([-1 -1 0; 1 -1 -1; 0 1 -1]; labels = ["plant", "herbivore", "predator"])
quasi_sign_stability(s)                 # (fraction, se)
press_sign_determinacy(s).sign          # sign of each response to a press
```

```@docs
SignPattern
quasi_sign_stability
press_sign_determinacy
```

## Applicability

```@docs
applicable_metrics
```
