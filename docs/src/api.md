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
linearise
PerturbationOptions
global_mortality_tolerance
local_mortality_tolerances
press_response
removal_responses
extinction_tolerance
```

## Applicability

```@docs
applicable_metrics
```
