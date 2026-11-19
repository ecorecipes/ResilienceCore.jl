# The 27 stability metrics of Domínguez-García, Dakos & Kéfi (2019), by name: 11 from the
# Jacobian at equilibrium and 16 from simulated perturbations of a community.

"""
    global_press_response(c; options = PerturbationOptions())

[`press_response`](@ref) to a press of `options.press` on every species' mortality at once, as
`(relative, l1, extinctions)`.
"""
function global_press_response(c; options::PerturbationOptions = PerturbationOptions())
    rel, l1, lost = press_response(c, fill(1 + options.press, length(state(c))); options = options)
    return (relative = rel, l1 = l1, extinctions = lost)
end

"""
    local_press_responses(c; options = PerturbationOptions())

[`press_response`](@ref) to a press of `options.press` on each present species' mortality in
turn, in the order of [`present_species`](@ref): a vector of `(relative, l1, extinctions)`.
"""
function local_press_responses(c; options::PerturbationOptions = PerturbationOptions())
    S = length(state(c))
    return map(present_species(c)) do i
        f = ones(S)
        f[i] = 1 + options.press
        rel, l1, lost = press_response(c, f; options = options)
        (relative = rel, l1 = l1, extinctions = lost)
    end
end

timekw(::EquilibriumJacobian{ContinuousTime}, dt) = (; dt = dt)
timekw(::EquilibriumJacobian{DiscreteTime}, dt) = (;)

"""
    jacobian_profile(op::EquilibriumJacobian; C = pulse_covariance(op), p = 2, dt = 0.1)

The 11 metrics of Domínguez-García et al. (2019) that depend only on the Jacobian, as a
`NamedTuple` in the paper's order:

| Field | Metric | Function |
|---|---|---|
| `R0` | reactivity | [`reactivity`](@ref) |
| `Is` | intrinsic stochastic invariability | [`stochastic_invariability`](@ref) |
| `Amax`, `tmax` | maximal amplification and its time | [`max_amplification`](@ref) |
| `Rinf` | asymptotic resilience | [`resilience`](@ref) |
| `MR0` | median reactivity | [`median_reactivity`](@ref) |
| `MAmax`, `Mtmax` | median maximal amplification and its time | [`median_max_amplification`](@ref) |
| `MRinf` | median resilience | [`median_resilience`](@ref) |
| `S` | community sensitivity `Σᵢ\\|Σⱼ sᵢⱼ\\|` | [`press_sensitivity`](@ref) |
| `mean_sij` | mean absolute sensitivity `Σᵢⱼ\\|sᵢⱼ\\|/N²` | [`total_sensitivity`](@ref) / N² |

Values keep their natural orientation (for example `R0 > 0` means the equilibrium amplifies
some perturbations); the paper flips signs so that larger means more stable. `C` is the pulse
covariance of the median metrics, `p` the norm of `R0` and `Amax`, and `dt` the time grid in
continuous time (ignored for maps).
"""
function jacobian_profile(op::EquilibriumJacobian; C::AbstractMatrix = pulse_covariance(op), p = 2, dt::Real = 0.1)
    kw = timekw(op, dt)
    N = size(op.J, 1)
    Amax, tmax = max_amplification(op; p = p, kw...)
    MAmax, Mtmax = median_max_amplification(op, C; kw...)
    return (R0 = reactivity(op; p = p), Is = stochastic_invariability(op), Amax = Amax, tmax = tmax,
            Rinf = resilience(op), MR0 = median_reactivity(op, C), MAmax = MAmax, Mtmax = Mtmax,
            MRinf = median_resilience(op, C; kw...), S = press_sensitivity(op), mean_sij = total_sensitivity(op) / N^2)
end

mean_of(x) = sum(x) / length(x)

"""
    stability_profile([rng], c; options = PerturbationOptions(), p = 2, dt = 0.1)

All 27 stability metrics of Domínguez-García et al. (2019) for a community that implements the
community interface ([`state`](@ref), [`present_species`](@ref), [`with_mortality`](@ref),
[`without_species`](@ref), [`jacobian`](@ref)), as a `NamedTuple` in the paper's order. The
first 11 are [`jacobian_profile`](@ref) of [`linearise`](@ref)`(c)`, with pulses proportional to
the abundances. The other 16 summarise the simulated experiments over the `N` species present:

| Field | Metric | From |
|---|---|---|
| `TM` | global mortality tolerance | [`global_mortality_tolerance`](@ref) |
| `mean_TM`, `min_TM` | mean and least local mortality tolerance | [`local_mortality_tolerances`](@ref) |
| `mean_CE`, `max_CE` | mean and largest number of secondary extinctions after one removal | [`removal_responses`](@ref) |
| `TE` | fraction of species removed before half are lost | [`extinction_tolerance`](@ref) / N |
| `RM` | \\|relative change in total abundance\\| under a global press | [`global_press_response`](@ref) |
| `mean_RM`, `max_RM` | \\|mean\\| and largest \\|relative change\\| under local presses | [`local_press_responses`](@ref) |
| `SM` | L1 change under a global press | [`global_press_response`](@ref) |
| `mean_SM`, `max_SM` | mean and largest L1 change under local presses | [`local_press_responses`](@ref) |
| `mean_RE`, `max_RE` | \\|mean\\| and largest \\|relative change\\| after one removal | [`removal_responses`](@ref) |
| `mean_SE`, `max_SE` | mean and largest L1 change after one removal | [`removal_responses`](@ref) |

Tolerances (`TM`, `mean_TM`, `min_TM`, `TE`) grow with stability; every other simulated metric
shrinks with it. The paper's `⟨CE⟩` and `CE_max` are `N/(N − 1)` times these counts, through
its normalisation; within a richness the ranks are the same. This runs every perturbation
experiment, so it is the expensive call; `rng` drives the random removal sequences of `TE`.
"""
function stability_profile(rng::AbstractRNG, c; options::PerturbationOptions = PerturbationOptions(), p = 2,
                           dt::Real = 0.1)
    N = richness(c)
    jac = jacobian_profile(linearise(c); C = pulse_covariance(abundances(c)), p = p, dt = dt)
    loc = local_mortality_tolerances(c; options = options)
    g = global_press_response(c; options = options)
    lp = local_press_responses(c; options = options)
    rem = removal_responses(c; options = options)
    secondary = [N - 1 - r.survivors for r in rem]
    lrel, ll1 = [x.relative for x in lp], [x.l1 for x in lp]
    rrel, rl1 = [r.relative for r in rem], [r.l1 for r in rem]
    return merge(jac,
                 (TM = global_mortality_tolerance(c; options = options), mean_TM = mean_of(loc),
                  min_TM = minimum(loc), mean_CE = mean_of(secondary),
                  TE = extinction_tolerance(rng, c; options = options) / N, max_CE = maximum(secondary),
                  RM = abs(g.relative), mean_RM = abs(mean_of(lrel)), max_RM = maximum(abs, lrel), SM = g.l1,
                  mean_SM = mean_of(ll1), max_SM = maximum(ll1), mean_RE = abs(mean_of(rrel)),
                  max_RE = maximum(abs, rrel), mean_SE = mean_of(rl1), max_SE = maximum(rl1)))
end
stability_profile(c; kwargs...) = stability_profile(Random.default_rng(), c; kwargs...)
