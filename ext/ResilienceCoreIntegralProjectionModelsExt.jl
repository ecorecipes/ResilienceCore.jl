module ResilienceCoreIntegralProjectionModelsExt

# linearise for integral projection models: the materialised kernel of a density-independent,
# deterministic IPM, with the mesh weights of its domain. Documented on the `linearise` binding
# in ResilienceCore.

using IntegralProjectionModels: IPMProblem, MegaKernel, materialize
using ResilienceCore: ResilienceCore, ProjectionMatrix, mesh_weights
using StructuredPopulationCore: AbstractStateDomain, DensityIndependent, Deterministic

# A MegaKernel stacks its states in the order of `keys(states)`, as `materialize` does.
kernel_weights(k::MegaKernel, domain) = reduce(vcat, [mesh_weights(k.states[s]) for s in keys(k.states)])
kernel_weights(k, domain::AbstractStateDomain) = mesh_weights(domain)
kernel_weights(k, domain) = throw(ArgumentError("cannot find the mesh of an IPMProblem with domain $(typeof(domain))"))

function ResilienceCore.linearise(prob::IPMProblem)
    prob.density isa DensityIndependent ||
        throw(ArgumentError("a density-dependent IPM has no single kernel; linearise its map at a fixed point " *
                            "instead, with `linearise(to_discrete_problem(prob), n_star)`"))
    prob.stochasticity isa Deterministic ||
        throw(ArgumentError("a stochastic IPM has no single kernel; use its stochastic growth rate"))
    K = Matrix{Float64}(materialize(prob.kernel))
    w = kernel_weights(prob.kernel, prob.domain)
    length(w) == size(K, 1) ||
        throw(ArgumentError("the kernel has $(size(K, 1)) states but the domain $(length(w)) mesh points"))
    return ProjectionMatrix(K; weights = w)
end

end
