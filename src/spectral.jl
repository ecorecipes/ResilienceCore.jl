# Asymptotic metrics: growth and return rates.

spectral_radius(M) = maximum(abs, eigvals(M))
spectral_abscissa(M) = maximum(real, eigvals(M))

"""
    growth_rate(op::ProjectionOperator)
    growth_rate(p::PeriodicOperator{<:ProjectionOperator})

The asymptotic growth rate of a density-independent model: the dominant eigenvalue λ of a
projection matrix (per step), the spectral abscissa of a generator (per unit time), or
`λ(period map)^(1/period)` for a periodic model. Not defined for equilibria, which do not grow.
"""
growth_rate(op::ProjectionOperator{DiscreteTime}) = spectral_radius(op.A)
growth_rate(op::ProjectionOperator{ContinuousTime}) = spectral_abscissa(op.A)
growth_rate(p::PeriodicOperator{<:ProjectionOperator}) = growth_rate(monodromy(p))^(1 / period(p))

"""
    resilience(op)

The asymptotic rate at which perturbations die out (larger is more resilient), per unit time or
per step:

- continuous-time equilibrium: `−Re λ₁(J)` (Pimm & Lawton);
- discrete-time fixed point: `−ln ρ(J)`;
- projection matrix: `ln(λ₁/|λ₂|)`, the log damping ratio, the rate of convergence to the
  stable structure (Caswell 2001);
- generator: `Re λ₁ − Re λ₂`;
- periodic model: the rate of its period map divided by the period;
- Markov kernel: `−ln |λ₂|`, the rate of convergence to the stationary distribution.

For an equilibrium, a negative value means it is unstable.
"""
resilience(op::EquilibriumJacobian{ContinuousTime}) = -spectral_abscissa(op.J)
resilience(op::EquilibriumJacobian{DiscreteTime}) = -log(spectral_radius(op.J))
function resilience(op::ProjectionOperator{DiscreteTime})
    m = sort(abs.(eigvals(op.A)); rev = true)
    return length(m) < 2 || m[2] == 0 ? Inf : log(m[1] / m[2])
end
function resilience(op::ProjectionOperator{ContinuousTime})
    r = sort(real.(eigvals(op.A)); rev = true)
    return length(r) < 2 ? Inf : r[1] - r[2]
end
resilience(p::PeriodicOperator) = resilience(monodromy(p)) / period(p)

"""
    isstable(op)

Whether an equilibrium or fixed point is locally asymptotically stable (`resilience(op) > 0`).
Not defined for density-independent models, whose populations grow or decline.
"""
isstable(op::EquilibriumJacobian) = resilience(op) > 0
isstable(p::PeriodicOperator{<:EquilibriumJacobian}) = resilience(p) > 0
