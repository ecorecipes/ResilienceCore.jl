# Qualitative (sign-pattern) stability: what can be said about a Jacobian when only the signs
# of its entries are known.

"""
    SignPattern(S; labels = nothing, unknown = falses(size(S)))

The sign pattern of a Jacobian: `S[i, j]` is the sign (−1, 0 or +1) of the effect of `j` on
`i`, as in a Jacobian or community matrix. `unknown[i, j]` marks a link that exists but whose
sign is not known; its entry in `S` is ignored. A zero diagonal means no self-regulation.
`labels` optionally names the variables (species, factors).

The metrics draw matrices with this pattern and magnitudes uniform on (0, 1). A link of
unknown sign is drawn as +1 or −1 with equal probability in each draw.
"""
struct SignPattern
    S::Matrix{Int8}
    unknown::BitMatrix
    labels::Union{Nothing,Vector{String}}
    function SignPattern(S::AbstractMatrix; labels = nothing, unknown::AbstractMatrix{Bool} = falses(size(S)))
        n = checksquare(S, "S")
        all(in((-1, 0, 1)), S) || throw(ArgumentError("a sign pattern has entries −1, 0 or +1"))
        size(unknown) == size(S) || throw(DimensionMismatch("unknown must have the size of S, $(size(S))"))
        labels === nothing || length(labels) == n ||
            throw(DimensionMismatch("$(length(labels)) labels for $n variables"))
        u = BitMatrix(unknown)
        s = Matrix{Int8}(S)
        s[u] .= 0
        return new(s, u, labels === nothing ? nothing : string.(collect(labels)))
    end
end

Base.size(s::SignPattern) = size(s.S)
Base.size(s::SignPattern, d) = size(s.S, d)

function Base.show(io::IO, s::SignPattern)
    print(io, "SignPattern(", size(s, 1), " variables, ", count(!iszero, s.S) + count(s.unknown), " links")
    any(s.unknown) && print(io, ", ", count(s.unknown), " of unknown sign")
    print(io, ")")
end

# A matrix with pattern `s` and magnitudes uniform on (0, 1). The magnitudes are drawn first,
# so a pattern without unknown links consumes the random stream exactly as `S .* rand(...)`.
function _draw(rng::AbstractRNG, s::SignPattern)
    M = rand(rng, size(s.S)...)
    J = s.S .* M
    if any(s.unknown)
        k = findall(s.unknown)
        J[k] .= M[k] .* ifelse.(rand(rng, Bool, length(k)), 1.0, -1.0)
    end
    return J
end

"""
    quasi_sign_stability([rng], s::SignPattern; samples = 10_000)

The fraction of matrices with sign pattern `s` and magnitudes uniform on (0, 1) that are
stable (every eigenvalue with negative real part): quasi-sign-stability (Allesina & Pascual
2008). Links of unknown sign take either sign with equal probability. Returns
`(fraction, se)`.
"""
function quasi_sign_stability(rng::AbstractRNG, s::SignPattern; samples::Integer = 10_000)
    stable = count(_ -> maximum(real, eigvals(_draw(rng, s))) < 0, 1:samples)
    p = stable / samples
    return (fraction = p, se = sqrt(p * (1 - p) / samples))
end
quasi_sign_stability(s::SignPattern; kwargs...) = quasi_sign_stability(Random.default_rng(), s; kwargs...)

"""
    press_sign_determinacy([rng], s::SignPattern; samples = 10_000)

Qualitative press-perturbation analysis (Dambacher et al. 2003; the simulation approach of
QPress): draw matrices with sign pattern `s` and magnitudes uniform on (0, 1), keep the stable
ones, and for each compute the equilibrium response to a sustained increase in each
variable's growth, `−J⁻¹`. Returns `(sign, determinacy, stable, labels)`: the most common
sign of each response (row: responding variable, column: pressed variable; 0 when the
response is structurally zero, below 1e-10 of the largest), the fraction of stable draws with
that outcome, the fraction of draws that were stable, and the pattern's labels.
"""
function press_sign_determinacy(rng::AbstractRNG, s::SignPattern; samples::Integer = 10_000)
    n = size(s, 1)
    pos, neg, kept = zeros(Int, n, n), zeros(Int, n, n), 0
    for _ in 1:samples
        J = _draw(rng, s)
        maximum(real, eigvals(J)) < 0 || continue
        R = -inv(J)
        tol = 1e-10 * maximum(abs, R)
        pos .+= R .> tol
        neg .+= R .< -tol
        kept += 1
    end
    kept == 0 && return (sign = zeros(Int8, n, n), determinacy = fill(NaN, n, n), stable = 0.0, labels = s.labels)
    zero_ = kept .- pos .- neg
    sgn = [zero_[i] ≥ max(pos[i], neg[i]) ? Int8(0) : Int8(sign(pos[i] - neg[i])) for i in CartesianIndices(pos)]
    det = max.(pos, neg, zero_) ./ kept
    return (sign = sgn, determinacy = det, stable = kept / samples, labels = s.labels)
end
press_sign_determinacy(s::SignPattern; kwargs...) = press_sign_determinacy(Random.default_rng(), s; kwargs...)
