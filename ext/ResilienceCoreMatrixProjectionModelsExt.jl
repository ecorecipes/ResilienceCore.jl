module ResilienceCoreMatrixProjectionModelsExt

# linearise for matrix projection models: a density-independent, deterministic MPM is a
# ProjectionMatrix. Documented on the `linearise` binding in ResilienceCore.

using MatrixProjectionModels: MatrixProjectionModel, MPMProblem
using ResilienceCore: ResilienceCore, ProjectionMatrix
using StructuredPopulationCore: DensityIndependent, Deterministic

ResilienceCore.linearise(m::MatrixProjectionModel) = ProjectionMatrix(m.A)

function ResilienceCore.linearise(prob::MPMProblem)
    prob.density isa DensityIndependent ||
        throw(ArgumentError("a density-dependent MPM has no single projection matrix; linearise its map at a " *
                            "fixed point instead, with `linearise(::DiscreteProblem, n_star)`"))
    prob.stochasticity isa Deterministic ||
        throw(ArgumentError("a stochastic MPM has no single projection matrix; use its stochastic growth rate"))
    M = prob.matrix
    M isa MatrixProjectionModel && return ProjectionMatrix(M.A)
    M isa AbstractMatrix && return ProjectionMatrix(M)
    throw(ArgumentError("cannot read a projection matrix from an MPMProblem holding a $(typeof(M))"))
end

end
