"""
    ResilienceCore

Resilience and stability metrics that dispatch on what a matrix *means*. The same numbers can be
the Jacobian of a flow at an equilibrium, the Jacobian of a map at a fixed point, a
density-independent projection matrix or generator, the period map of a seasonal model, or a
Markov kernel. Each metric has methods only for the objects on which it is defined, so asking
for, say, the growth rate of an equilibrium Jacobian is a `MethodError` rather than a number
without meaning. [`applicable_metrics`](@ref) lists what can be computed for an object.

Time semantics reuse `DiscreteTime` and `ContinuousTime` from StructuredPopulationCore, the
shared core of the ecorecipes projection-model packages.
"""
module ResilienceCore

using LinearAlgebra: I, Diagonal, Symmetric, eigen, eigmax, eigvals, inv, issuccess, kron, lu, norm, opnorm, svdvals,
                     sylvester, tr
using Random: Random, AbstractRNG
using StructuredPopulationCore: AbstractTimeSemantics, ContinuousTime, DiscreteTime, ContinuousDomain, DiscreteDomain,
                                step_size, n_states

export AbstractTimeSemantics, ContinuousTime, DiscreteTime
export LinearisedDynamics, EquilibriumJacobian, ContinuousJacobian, DiscreteJacobian, ProjectionOperator,
       ProjectionMatrix, Generator, PeriodicOperator, MarkovKernel
export timesemantics, operator, monodromy, period
export growth_rate, resilience, isstable, reactivity, amplification_envelope, max_amplification, propagator
export pulse_covariance, median_reactivity, median_displacement, median_max_amplification, median_resilience
export sensitivity_matrix, press_sensitivity, total_sensitivity, stochastic_invariability, growth_sensitivity,
       elasticity
export stationary_distribution, mixing_time
export applicable_metrics
export state, present_species, with_mortality, without_species, jacobian, richness, abundances, linearise
export mesh_weights, affine_equilibrium, periodic_generator
export PerturbationOptions, global_mortality_tolerance, local_mortality_tolerances, press_response, removal_responses,
       extinction_tolerance

include("operators.jl")
include("spectral.jl")
include("transient.jl")
include("median.jl")
include("press.jl")
include("markov.jl")
include("applicable.jl")
include("community.jl")
include("linearise.jl")

end
