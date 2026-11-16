# Simulation-based resilience: experiments on any community model that can be perturbed and
# brought back to a steady state. A model package implements the interface below for its own
# community type; the experiments then work unchanged.

"""
    state(c)

The abundance of every species the community started with, in a fixed order, zero for those
absent. Part of the community interface.
"""
function state end

"""
    present_species(c)

Indices (into [`state`](@ref)) of the species present. Part of the community interface.
"""
function present_species end

"""
    with_mortality(c, factors; newton = true, steady = nothing)

The community after setting every species' mortality to `factors` times its baseline (a vector
over all species, in the order of [`state`](@ref)) and moving it to the new steady state from
its current state. `newton` allows a model to try a Newton step first. `steady` passes
model-specific steady-state settings (`nothing`: the model's default for perturbations). Part
of the community interface.
"""
function with_mortality end

"""
    without_species(c, i; steady = nothing)

The community after removing species `i` and moving to the new steady state. Part of the
community interface.
"""
function without_species end

"""
    jacobian(c)

The Jacobian of a community model at its current state, over the species present. Model
packages add methods for their own types.
"""
function jacobian end

"""
    richness(c)

The number of species present.
"""
richness(c) = length(present_species(c))

"""
    abundances(c)

The abundances of the species present, `state(c)[present_species(c)]`.
"""
abundances(c) = state(c)[present_species(c)]

"""
    linearise(c)

The `ContinuousJacobian` (an [`EquilibriumJacobian`](@ref)) of a continuous-time community at
its current state, with the abundances of the species present as its equilibrium.
"""
linearise(c) = ContinuousJacobian(jacobian(c), abundances(c))

"""
    PerturbationOptions(; press = 0.1, global_start = 1.001, step = 0.1, coarse = 10,
                          max_steps = 2000, sequences = 10, steady = nothing)

Settings of the simulated perturbation experiments:

- `press`: the relative size of the mortality press for resistance and sensitivity (10%).
- `global_start`, `step`: the first global mortality multiplier and the grid step of the
  tolerance searches.
- `coarse`: the first jump of the tolerance search, in grid steps. It doubles while no species
  dies and halves when one does, down to single grid steps.
- `max_steps`: the largest number of grid steps tried.
- `sequences`: the number of random removal sequences for [`extinction_tolerance`](@ref).
- `steady`: model-specific steady-state settings, passed through (`nothing` uses the model's
  own default for perturbations).

The defaults follow Domínguez-García et al. (2019).
"""
Base.@kwdef struct PerturbationOptions
    press::Float64 = 0.1
    global_start::Float64 = 1.001
    step::Float64 = 0.1
    coarse::Int = 10
    max_steps::Int = 2000
    sequences::Int = 10
    steady::Any = nothing
end

# Smallest grid index k ≥ k0 at which mortality multipliers `factors(k)` drive some species
# extinct. Each level is reached from the last extinction-free state. Jumps start at
# `opts.coarse` grid steps and double while no species dies; when a jump kills one it is
# halved and retried, so the deciding step is always a single grid increment. This assumes a
# jump that causes no extinction would cause none if taken in smaller steps.
function first_extinction(c, factors::Function, k0::Integer, opts::PerturbationOptions)
    n0 = richness(c)
    attempt(base, k) = (t = with_mortality(base, factors(k); steady = opts.steady); (richness(t) == n0, t))
    ok, cur = attempt(c, k0)
    ok || return k0
    lo, jump = k0, max(opts.coarse, 1)
    while lo < opts.max_steps
        hi = min(lo + jump, opts.max_steps)
        ok, trial = attempt(cur, hi)
        if ok
            cur, lo = trial, hi
            jump *= 2
        elseif hi == lo + 1
            return hi
        else
            jump = max((hi - lo) ÷ 2, 1)
        end
    end
    return nothing
end

"""
    global_mortality_tolerance(c; options = PerturbationOptions())

The mortality multiplier `global_start + step k` (k = 0, 1, …), applied to every species, at
which the first extinction occurs (Wootton & Stouffer 2016). `NaN` if none occurs within
`max_steps`.
"""
function global_mortality_tolerance(c; options::PerturbationOptions = PerturbationOptions())
    S = length(state(c))
    f(k) = options.global_start + options.step * k
    k = first_extinction(c, k -> fill(f(k), S), 0, options)
    return k === nothing ? NaN : f(k)
end

"""
    local_mortality_tolerances(c; options = PerturbationOptions())

For each species present, the smallest relative increase `step k` (k = 1, 2, …) in its own
mortality that causes any extinction, in the order of [`present_species`](@ref).
"""
function local_mortality_tolerances(c; options::PerturbationOptions = PerturbationOptions())
    S = length(state(c))
    out = Float64[]
    for i in present_species(c)
        factors(k) = (f = ones(S); f[i] = 1 + options.step * k; f)
        k = first_extinction(c, factors, 1, options)
        push!(out, k === nothing ? NaN : options.step * k)
    end
    return out
end

"""
    press_response(c, factors; options = PerturbationOptions(), newton = true)

Set mortality multipliers `factors` (a vector over all species), move to the new steady state,
and return `(relative_change, l1_change, extinctions)`. These are `(B₀ − B′)/B₀` for total
abundance, the L1 norm of the species-wise change divided by `B₀`, and the number of species
lost.
"""
function press_response(c, factors::AbstractVector; options::PerturbationOptions = PerturbationOptions(),
                        newton::Bool = true)
    p = with_mortality(c, factors; newton = newton, steady = options.steady)
    x0, x = state(c), state(p)
    B0 = sum(x0)
    return ((B0 - sum(x)) / B0, sum(abs, x0 .- x) / B0, richness(c) - richness(p))
end

"""
    removal_responses(c; options = PerturbationOptions())

Remove each species present in turn and move to the new steady state. Returns a vector of
`(survivors, relative, l1)`: the number of species left (at most N − 1) and the relative and
L1 changes in abundance, as in [`press_response`](@ref).
"""
function removal_responses(c; options::PerturbationOptions = PerturbationOptions())
    out = NamedTuple{(:survivors, :relative, :l1),Tuple{Int,Float64,Float64}}[]
    x0 = state(c)
    B0 = sum(x0)
    for i in present_species(c)
        r = without_species(c, i; steady = options.steady)
        x = state(r)
        push!(out, (survivors = richness(r), relative = (B0 - sum(x)) / B0, l1 = sum(abs, x0 .- x) / B0))
    end
    return out
end

"""
    extinction_tolerance([rng], c; options = PerturbationOptions())

Robustness (Dunne et al. 2002): the number of random primary removals, each followed by a
return to steady state, needed before at most half of the original species remain, averaged
over `options.sequences` random sequences.
"""
function extinction_tolerance(rng::AbstractRNG, c; options::PerturbationOptions = PerturbationOptions())
    n0 = richness(c)
    total = 0
    for _ in 1:options.sequences
        cur, removed = c, 0
        while richness(cur) > n0 / 2
            cur = without_species(cur, rand(rng, present_species(cur)); steady = options.steady)
            removed += 1
        end
        total += removed
    end
    return total / options.sequences
end
extinction_tolerance(c; kwargs...) = extinction_tolerance(Random.default_rng(), c; kwargs...)
