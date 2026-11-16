# Operators that carry their meaning in their type.

"""
    LinearisedDynamics{T<:AbstractTimeSemantics}

A matrix together with what it means for the dynamics, in continuous (`ContinuousTime`) or
discrete (`DiscreteTime`) time. Subtypes:

- [`EquilibriumJacobian`](@ref): perturbations about an equilibrium or fixed point, which decay
  to zero when it is stable;
- [`ProjectionOperator`](@ref): density-independent growth, where perturbations of the
  population *structure* decay relative to the dominant mode.
"""
abstract type LinearisedDynamics{T<:AbstractTimeSemantics} end

"""
    timesemantics(op)

`ContinuousTime()` or `DiscreteTime()`.
"""
timesemantics(::LinearisedDynamics{T}) where {T} = T()

function checksquare(M, name)
    size(M, 1) == size(M, 2) || throw(ArgumentError("$name must be square, got size $(size(M))"))
    return size(M, 1)
end

"""
    EquilibriumJacobian{T}(J, equilibrium = Float64[])
    ContinuousJacobian(J, equilibrium = Float64[])
    DiscreteJacobian(J, equilibrium = Float64[])

The Jacobian `J` of a non-linear model at an equilibrium (`ContinuousTime`, `dx/dt = f(x)`) or a
fixed point (`DiscreteTime`, `x′ = f(x)`). Small perturbations evolve as `exp(J t)` or `Jᵗ` and
decay to zero when the point is stable. `equilibrium` (the state at the point, e.g. species
abundances) is needed only for the median metrics, whose default pulses are proportional to it.
"""
struct EquilibriumJacobian{T<:AbstractTimeSemantics} <: LinearisedDynamics{T}
    J::Matrix{Float64}
    equilibrium::Vector{Float64}
    function EquilibriumJacobian{T}(J::AbstractMatrix, equilibrium::AbstractVector = Float64[]) where {T<:AbstractTimeSemantics}
        n = checksquare(J, "J")
        isempty(equilibrium) || length(equilibrium) == n ||
            throw(ArgumentError("equilibrium has $(length(equilibrium)) entries for a $n × $n Jacobian"))
        return new{T}(Matrix{Float64}(J), Vector{Float64}(equilibrium))
    end
end
const ContinuousJacobian = EquilibriumJacobian{ContinuousTime}
const DiscreteJacobian = EquilibriumJacobian{DiscreteTime}

"""
    ProjectionOperator{T}(A; weights = ones(n))
    ProjectionMatrix(A; weights)    # DiscreteTime: n(t + 1) = A n(t)
    Generator(Q; weights)           # ContinuousTime: dn/dt = Q n

A density-independent linear model: a non-negative projection matrix (matrix projection model,
or a discretised integral projection kernel) or a generator with non-negative off-diagonal
entries (finite-state continuous-time model). There is no equilibrium; the population grows at
the rate [`growth_rate`](@ref), and the stability metrics describe how perturbations of its
structure decay *relative to that growth*.

`weights` are the mesh weights of a discretised continuous state (an IPM's bin widths), so that
norms approximate integrals and the metrics do not depend on the mesh. Use ones for a matrix
model.
"""
struct ProjectionOperator{T<:AbstractTimeSemantics} <: LinearisedDynamics{T}
    A::Matrix{Float64}
    weights::Vector{Float64}
    function ProjectionOperator{T}(A::AbstractMatrix; weights::AbstractVector = ones(size(A, 1))) where {T<:AbstractTimeSemantics}
        n = checksquare(A, "A")
        length(weights) == n && all(>(0), weights) ||
            throw(ArgumentError("weights must be $n positive numbers"))
        if T === DiscreteTime
            all(≥(0), A) || throw(ArgumentError("a projection matrix must be non-negative"))
        else
            all(A[i, j] ≥ 0 for i in 1:n, j in 1:n if i != j) ||
                throw(ArgumentError("a generator must have non-negative off-diagonal entries"))
        end
        return new{T}(Matrix{Float64}(A), Vector{Float64}(weights))
    end
end
const ProjectionMatrix = ProjectionOperator{DiscreteTime}
const Generator = ProjectionOperator{ContinuousTime}

"""
    operator(op)

The matrix of an operator.
"""
operator(op::EquilibriumJacobian) = op.J
operator(op::ProjectionOperator) = op.A

weights(op::EquilibriumJacobian) = ones(size(op.J, 1))
weights(op::ProjectionOperator) = op.weights

# the same kind of operator with another matrix
rebuild(op::EquilibriumJacobian{T}, M) where {T} = EquilibriumJacobian{T}(M, op.equilibrium)
rebuild(op::ProjectionOperator{T}, M) where {T} = ProjectionOperator{T}(M; weights = op.weights)

"""
    PeriodicOperator(parts)

A periodic (e.g. seasonal) discrete-time model: the operators `parts[1], parts[2], …` are
applied in turn, all of the same kind. Its metrics are those of the period map
[`monodromy`](@ref), expressed per step.
"""
struct PeriodicOperator{K<:LinearisedDynamics{DiscreteTime}}
    parts::Vector{K}
    function PeriodicOperator(parts::AbstractVector{K}) where {K<:LinearisedDynamics{DiscreteTime}}
        isempty(parts) && throw(ArgumentError("a periodic operator needs at least one part"))
        isconcretetype(K) || throw(ArgumentError("all parts must be operators of the same kind, got $K"))
        n = size(operator(first(parts)), 1)
        all(p -> size(operator(p), 1) == n, parts) || throw(ArgumentError("all parts must have the same size"))
        return new{K}(collect(parts))
    end
end

"""
    period(p::PeriodicOperator)

The number of steps in a period.
"""
period(p::PeriodicOperator) = length(p.parts)

"""
    monodromy(p::PeriodicOperator)

The period map `A_p ⋯ A₂ A₁`, as an operator of the same kind as the parts.
"""
monodromy(p::PeriodicOperator) =
    rebuild(first(p.parts), foldl((M, part) -> operator(part) * M, p.parts[2:end]; init = operator(first(p.parts))))

"""
    MarkovKernel(P)

A row-stochastic transition matrix (`P[i, j]` is the probability of moving from `i` to `j`),
which moves *distributions*, not abundances: its dominant eigenvalue is always 1, so
growth-based and equilibrium-based metrics do not apply. See [`resilience`](@ref),
[`mixing_time`](@ref).
"""
struct MarkovKernel
    P::Matrix{Float64}
    function MarkovKernel(P::AbstractMatrix)
        checksquare(P, "P")
        all(≥(0), P) && all(r -> isapprox(sum(r), 1; atol = 1e-10), eachrow(P)) ||
            throw(ArgumentError("P must be non-negative with rows summing to 1"))
        return new(Matrix{Float64}(P))
    end
end
