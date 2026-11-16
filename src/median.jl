# Median metrics over random pulses (Arnoldi et al. 2016, 2018), for equilibria and fixed points.

"""
    pulse_covariance(B)
    pulse_covariance(op::EquilibriumJacobian)

The covariance `Diagonal(B.^2)` of pulses proportional to the equilibrium state `B`
(Arnoldi et al. 2018; Domínguez-García et al. 2019).
"""
pulse_covariance(B::AbstractVector) = Diagonal(B .^ 2)
function pulse_covariance(op::EquilibriumJacobian)
    isempty(op.equilibrium) &&
        throw(ArgumentError("median metrics need the equilibrium of the Jacobian, or an explicit pulse covariance C"))
    return pulse_covariance(op.equilibrium)
end

"""
    median_displacement(op, C, t)

The median size of a pulse of covariance `C` after time `t` relative to its initial size,
`sqrt(Tr(C Φᵀ Φ) / Tr C)` with `Φ` the propagator (`exp(J t)`, or `Jᵗ` for integer `t`).
"""
function median_displacement(op::EquilibriumJacobian, C::AbstractMatrix, t)
    Φ = propagator(op, t; p = 2)
    return sqrt(max(tr(C * (Φ' * Φ)), 0.0) / tr(C))
end

"""
    median_reactivity(op, C = pulse_covariance(op))

The median initial rate of displacement: `Tr(C J) / Tr C` in continuous time, and
`ln median_displacement(op, C, 1)` per step in discrete time (the two agree when `J_d = exp(J Δt)`
and `Δt → 0`).
"""
median_reactivity(op::EquilibriumJacobian{ContinuousTime}, C::AbstractMatrix = pulse_covariance(op)) =
    tr(C * op.J) / tr(C)
median_reactivity(op::EquilibriumJacobian{DiscreteTime}, C::AbstractMatrix = pulse_covariance(op)) =
    log(median_displacement(op, C, 1))

"""
    median_max_amplification(op, C = pulse_covariance(op); dt = 0.1, threshold = 0.005, tmax)

`(MAmax, Mtmax)`: the maximum of [`median_displacement`](@ref) and when it is reached, scanning
until the displacement falls below `threshold` (time grid `dt` in continuous time, integer
steps in discrete time). The scan does not stop when the displacement first falls: a median
pulse can shrink and then grow.
"""
function median_max_amplification(op::EquilibriumJacobian{ContinuousTime}, C::AbstractMatrix = pulse_covariance(op);
                                  dt::Real = 0.1, threshold::Real = 0.005, tmax::Real = 1e4)
    return median_scan(C, exp(op.J * dt), dt, floor(Int, tmax / dt), threshold)
end
function median_max_amplification(op::EquilibriumJacobian{DiscreteTime}, C::AbstractMatrix = pulse_covariance(op);
                                  threshold::Real = 0.005, tmax::Integer = 10^4)
    best, k = median_scan(C, op.J, 1, tmax, threshold)
    return (best, round(Int, k))
end
function median_scan(C, step, dt, nmax, threshold)
    Φ, trC = Matrix(step), tr(C)
    best, tbest = 1.0, 0.0
    for k in 1:nmax
        M = sqrt(max(tr(C * (Φ' * Φ)), 0.0) / trC)
        M > best && ((best, tbest) = (M, k * dt))
        M < threshold && break
        Φ = Φ * step
    end
    return (best, round(tbest; digits = 10))
end

"""
    median_resilience(op, C = pulse_covariance(op); threshold = 0.005, dt = 0.1)

The median return rate `−ln(median_displacement²) / (2t)`, at the first time on the grid (step
`dt`, or integer steps) at which the displacement is below `threshold`. `NaN` if the
equilibrium is unstable.
"""
function median_resilience(op::EquilibriumJacobian{ContinuousTime}, C::AbstractMatrix = pulse_covariance(op);
                           threshold::Real = 0.005, dt::Real = 0.1)
    isstable(op) || return NaN
    return first_below(k -> median_displacement(op, C, k * dt), threshold, dt)
end
function median_resilience(op::EquilibriumJacobian{DiscreteTime}, C::AbstractMatrix = pulse_covariance(op);
                           threshold::Real = 0.005)
    isstable(op) || return NaN
    return first_below(k -> median_displacement(op, C, k), threshold, 1)
end
# doubling then bisection on the grid index for the first displacement below the threshold
function first_below(disp, threshold, dt)
    hi = 1
    while disp(hi) ≥ threshold
        hi *= 2
        hi > 2^40 && return NaN
    end
    lo = hi ÷ 2
    while hi - lo > 1
        mid = (lo + hi) ÷ 2
        disp(mid) < threshold ? (hi = mid) : (lo = mid)
    end
    return -log(disp(hi)^2) / (2 * hi * dt)
end
