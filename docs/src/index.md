# ResilienceCore.jl

```@docs
ResilienceCore
```

## Which metric for which object

| | `ContinuousJacobian` | `DiscreteJacobian` | `ProjectionMatrix` | `Generator` | `PeriodicOperator` | `MarkovKernel` |
|---|---|---|---|---|---|---|
| `growth_rate` | — | — | λ = ρ(A) | spectral abscissa | λ(period map)^(1/p) (projections) | — |
| `resilience` | −Re λ₁ | −ln ρ(J) | ln(λ₁/\|λ₂\|) | Re λ₁ − Re λ₂ | rate of the period map / p | −ln \|λ₂\| |
| `isstable` | ✓ | ✓ | — | — | ✓ (Jacobians) | — |
| `reactivity` | μ_p(J) | ln ‖J‖_p | ln ‖A/λ‖_p | μ_p(Q − sI) | — | — |
| `max_amplification` | max ‖e^{Jt}‖ | max ‖Jᵗ‖ | max ‖(A/λ)ᵗ‖ | max ‖e^{(Q−sI)t}‖ | — | — |
| median metrics | ✓ | ✓ | — | — | — | — |
| `press_sensitivity`, `sensitivity_matrix` | −J⁻¹ | (I − J)⁻¹ | — | — | — | — |
| `stochastic_invariability` | Lyapunov | Stein | — | — | — | — |
| `growth_sensitivity`, `elasticity` | — | — | ✓ | sensitivity | — | — |
| `stationary_distribution`, `mixing_time` | — | — | — | — | — | ✓ |

A dash is a `MethodError`; [`applicable_metrics`](@ref) gives the list for any object.
Norms default to the Euclidean norm for Jacobians (Neubert & Caswell 1997). For projections the
default is the 1-norm, i.e. total population (Stott et al. 2011), weighted by the mesh for an
integral projection model.

The metrics are tied together by identities that the tests check:
- the metrics of `exp(J Δt)` per `Δt` equal those of `J`;
- a rescaled projection matrix has the same transients;
- growth rate and damping ratio agree with StructuredPopulationCore;
- a weighted integral-projection kernel gives the same metrics on any mesh.
