# Responses to sustained (press) and stochastic perturbations.

"""
    sensitivity_matrix(op::EquilibriumJacobian)

The change in the equilibrium per unit sustained press on each variable (Bender et al. 1984):
`−J⁻¹` for a flow and `(I − J)⁻¹` for a map. The two agree for `J_d = I + Δt J` (with the press
scaled by `Δt`).
"""
sensitivity_matrix(op::EquilibriumJacobian{ContinuousTime}) = -inv(op.J)
sensitivity_matrix(op::EquilibriumJacobian{DiscreteTime}) = inv(I - op.J)

"""
    press_sensitivity(op)

Community sensitivity `Σᵢ |Σⱼ sᵢⱼ|` of the [`sensitivity_matrix`](@ref) (Domínguez-García et al.
2019, SI eq. 15).
"""
press_sensitivity(op::EquilibriumJacobian) = sum(abs, sum(sensitivity_matrix(op); dims = 2))

"""
    total_sensitivity(op)

`Σᵢⱼ |sᵢⱼ|`, the summed absolute entries of the [`sensitivity_matrix`](@ref).
"""
total_sensitivity(op::EquilibriumJacobian) = sum(abs, sensitivity_matrix(op))

"""
    stochastic_invariability(op::EquilibriumJacobian; dense_max)

Intrinsic stochastic invariability (Arnoldi et al. 2016): the inverse of the largest stationary
variance response to white noise.

- **Flow:** `1 / (2‖Ĵ⁻¹‖₂)` with `Ĵ = J ⊗ I + I ⊗ J`, from the Lyapunov equation. For
  `n > dense_max` it is found by power iteration on Sylvester solves, in O(n²) memory.
- **Map:** `σ_min(I − J ⊗ J) / 2` per step, from the discrete Lyapunov (Stein) equation, which
  equals the flow value times `Δt` for `J_d = exp(J Δt)`, `Δt → 0`. Dense only; `n` must be at
  most `dense_max`.
"""
function stochastic_invariability(op::EquilibriumJacobian{ContinuousTime}; dense_max::Integer = 30, tol::Real = 1e-10,
                                  maxiter::Integer = 2000)
    J = op.J
    n = size(J, 1)
    if n ≤ dense_max
        In = Matrix{Float64}(I, n, n)
        return minimum(svdvals(kron(J, In) + kron(In, J))) / 2
    end
    Jt = Matrix(J')
    X = ones(n, n) ./ n
    λ = 0.0
    for _ in 1:maxiter
        Y = sylvester(J, Jt, -X)        # J Y + Y Jᵀ = X
        Z = sylvester(Jt, J, -Y)        # Jᵀ Z + Z J = Y
        λnew = norm(Z)
        X = Z ./ λnew
        abs(λnew - λ) ≤ tol * λnew && (λ = λnew; break)
        λ = λnew
    end
    return 1 / (2sqrt(λ))
end
function stochastic_invariability(op::EquilibriumJacobian{DiscreteTime}; dense_max::Integer = 60)
    J = op.J
    n = size(J, 1)
    n ≤ dense_max || throw(ArgumentError("discrete-time invariability is computed densely; n = $n > dense_max = $dense_max"))
    return minimum(svdvals(I - kron(J, J))) / 2
end

# dominant right (w) and left (v) eigenvectors
function dominant_vectors(A)
    F, G = eigen(A), eigen(Matrix(A'))
    w = real.(F.vectors[:, argmax(real.(F.values))])
    v = real.(G.vectors[:, argmax(real.(G.values))])
    return w, v
end

"""
    growth_sensitivity(op::ProjectionOperator)

The sensitivity of the growth rate to each matrix entry, `∂λ/∂aᵢⱼ = vᵢ wⱼ / ⟨v, w⟩`, with `w` the
stable structure and `v` the reproductive values (Caswell 2001). For a generator, the same for
the spectral abscissa.
"""
function growth_sensitivity(op::ProjectionOperator)
    w, v = dominant_vectors(op.A)
    return v * w' ./ (v' * w)
end

"""
    elasticity(op::ProjectionMatrix)

Proportional sensitivities `eᵢⱼ = (aᵢⱼ/λ) ∂λ/∂aᵢⱼ`; they sum to one.
"""
elasticity(op::ProjectionOperator{DiscreteTime}) = op.A .* growth_sensitivity(op) ./ growth_rate(op)
