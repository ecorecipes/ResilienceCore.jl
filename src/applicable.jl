# Which metrics are defined for which objects.

const METRICS = (growth_rate, resilience, isstable, reactivity, max_amplification, stochastic_invariability,
                 sensitivity_matrix, press_sensitivity, total_sensitivity, growth_sensitivity, elasticity,
                 median_reactivity, median_max_amplification, median_resilience, stationary_distribution,
                 mixing_time)

"""
    applicable_metrics(x)
    applicable_metrics(T::Type)

The metric functions defined for `x` (or for objects of type `T`). A metric that is not listed
throws a `MethodError` whose message names the ones that are.
"""
applicable_metrics(x) = applicable_metrics(typeof(x))
applicable_metrics(T::Type) = [f for f in METRICS if hasmethod(f, Tuple{T})]

const OPERATORS = Union{LinearisedDynamics,PeriodicOperator,MarkovKernel}

function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
        (exc.f in METRICS && !isempty(argtypes) && argtypes[1] <: OPERATORS) || return
        ok = applicable_metrics(argtypes[1])
        print(io, "\n`", nameof(exc.f), "` is not defined for a ", argtypes[1], ". Metrics defined for it: ",
              isempty(ok) ? "none" : join(nameof.(ok), ", "), ".")
    end
end
