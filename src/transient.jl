# Transient metrics: reactivity and amplification, on the operator standardised so that its
# asymptotic behaviour is neutral (equilibria as they are; projections relative to their growth).

default_norm(::EquilibriumJacobian) = 2   # Euclidean, as Neubert & Caswell (1997)
default_norm(::ProjectionOperator) = 1    # population size, as Stott et al. (2011)

# The matrix whose powers or exponential describe transients: J for an equilibrium; A/λ for a
# projection matrix; Q − s I for a generator.
base(op::EquilibriumJacobian) = op.J
base(op::ProjectionOperator{DiscreteTime}) = op.A ./ growth_rate(op)
base(op::ProjectionOperator{ContinuousTime}) = op.A - growth_rate(op) * I

function norm_weights(w, p)
    p == 1 && return w
    p == 2 && return sqrt.(w)
    p == Inf && return one.(w)
    throw(ArgumentError("p must be 1, 2 or Inf, got $p"))
end

# `base` in coordinates where the weighted p-norm is the plain p-norm
function standardised(op::LinearisedDynamics, p)
    d = norm_weights(weights(op), p)
    return Diagonal(d) * base(op) * Diagonal(1 ./ d)
end

# logarithmic norm μ_p(M), the initial growth rate of ‖exp(M t)‖_p
function lognorm(M, p)
    p == 2 && return eigmax(Symmetric((M + M') / 2))
    n = size(M, 1)
    p == 1 && return maximum(M[j, j] + sum(abs(M[i, j]) for i in 1:n if i != j; init = 0.0) for j in 1:n)
    p == Inf && return maximum(M[i, i] + sum(abs(M[i, j]) for j in 1:n if j != i; init = 0.0) for i in 1:n)
    throw(ArgumentError("p must be 1, 2 or Inf, got $p"))
end

"""
    reactivity(op; p = 2 for equilibria, 1 for projections)

The largest initial rate at which a perturbation can grow, per unit time or per step:

- continuous time: the logarithmic norm `μ_p` of `J` (for p = 2, the largest eigenvalue of
  `(J + Jᵀ)/2`; Neubert & Caswell 1997), or of `Q − s I` for a generator;
- discrete time: `ln ‖J‖_p`, or `ln ‖A/λ‖_p` for a projection matrix (Stott et al. 2011 use
  the 1-norm, the total population, without the logarithm).

Positive means some perturbation grows at first (the system is *reactive*).
"""
reactivity(op::LinearisedDynamics{ContinuousTime}; p = default_norm(op)) = lognorm(standardised(op, p), p)
reactivity(op::LinearisedDynamics{DiscreteTime}; p = default_norm(op)) = log(opnorm(standardised(op, p), p))

"""
    propagator(op, t; p)

The standardised propagator: `exp(M t)` in continuous time, `Mᵗ` (integer `t`) in discrete time,
with `M` the operator relative to its asymptotic growth, in weighted coordinates.
"""
propagator(op::LinearisedDynamics{ContinuousTime}, t::Real; p = default_norm(op)) = exp(standardised(op, p) * t)
propagator(op::LinearisedDynamics{DiscreteTime}, t::Integer; p = default_norm(op)) = standardised(op, p)^t

"""
    amplification_envelope(op, t; p)

`‖propagator(op, t)‖_p`: the largest factor by which any perturbation is amplified at time `t`
(Neubert & Caswell 1997), relative to asymptotic growth for a projection.
"""
amplification_envelope(op::LinearisedDynamics, t; p = default_norm(op)) = opnorm(propagator(op, t; p = p), p)

# the limit of the standardised propagator: 0 for a stable equilibrium, the projector onto the
# dominant mode for a projection
asymptote(op::EquilibriumJacobian, M) = zero(M)
function asymptote(op::ProjectionOperator, M)
    F, G = eigen(M), eigen(Matrix(M'))
    r = real.(F.vectors[:, argmax(real.(F.values))])
    l = real.(G.vectors[:, argmax(real.(G.values))])
    return r * l' ./ (l' * r)
end

"""
    max_amplification(op; p, dt = 0.1, tmax = 1e4, settle = 1e-8)

`(Amax, tmax)`: the maximum of [`amplification_envelope`](@ref) and when it is reached, on a time
grid of step `dt` (continuous time) or over integer steps (discrete time). If the operator is
not reactive in the chosen norm, the envelope never exceeds its initial value of 1 and the
result is `(1.0, 0)`. Otherwise the scan stops when the propagator is within `settle` of its
limit (zero for an equilibrium; the dominant-mode projector for a projection, whose norm is the
asymptotic amplification) or at `tmax`.

The amplification of an unstable equilibrium is unbounded, so the result is `(Inf, Inf)`. The
same is returned if the propagator overflows during the scan.
"""
function max_amplification(op::LinearisedDynamics{ContinuousTime}; p = default_norm(op), dt::Real = 0.1,
                           tmax::Real = 1e4, settle::Real = 1e-8)
    op isa EquilibriumJacobian && !isstable(op) && return (Inf, Inf)
    reactivity(op; p = p) ≤ 0 && return (1.0, 0.0)
    M = standardised(op, p)
    return scan(op, M, exp(M * dt), dt, floor(Int, tmax / dt), p, settle)
end
function max_amplification(op::LinearisedDynamics{DiscreteTime}; p = default_norm(op), tmax::Integer = 10^4,
                           settle::Real = 1e-8)
    op isa EquilibriumJacobian && !isstable(op) && return (Inf, Inf)
    reactivity(op; p = p) ≤ 0 && return (1.0, 0)
    M = standardised(op, p)
    best, k = scan(op, M, M, 1, tmax, p, settle)
    return isfinite(best) ? (best, round(Int, k)) : (Inf, Inf)
end

function scan(op, M, step, dt, nmax, p, settle)
    limit = asymptote(op, M)
    Φ = Matrix(step)
    best, tbest = 1.0, 0.0
    for k in 1:nmax
        all(isfinite, Φ) || return (Inf, Inf)
        A = opnorm(Φ, p)
        A > best && ((best, tbest) = (A, k * dt))
        opnorm(Φ - limit, p) < settle * max(1.0, opnorm(limit, p)) && break
        op isa EquilibriumJacobian && A < 1 && tbest > 0 && break
        Φ = Φ * step
    end
    return (best, round(tbest; digits = 10))
end
