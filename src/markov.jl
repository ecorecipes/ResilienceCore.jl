# Markov kernels: convergence of distributions.

function resilience(m::MarkovKernel)
    a = sort(abs.(eigvals(m.P)); rev = true)
    return length(a) < 2 || a[2] == 0 ? Inf : -log(a[2])
end

"""
    stationary_distribution(m::MarkovKernel)

The left eigenvector of `P` for eigenvalue 1, normalised to sum to one.
"""
function stationary_distribution(m::MarkovKernel)
    F = eigen(Matrix(m.P'))
    v = real.(F.vectors[:, argmin(abs.(F.values .- 1))])
    return v ./ sum(v)
end

"""
    mixing_time(m::MarkovKernel; ε = 0.25, tmax = 10^6)

The first step `t` at which the distribution from every starting state is within total
variation `ε` of the stationary distribution, or `nothing` if not reached by `tmax`.
"""
function mixing_time(m::MarkovKernel; ε::Real = 0.25, tmax::Integer = 10^6)
    stat = stationary_distribution(m)
    Pt = copy(m.P)
    for t in 1:tmax
        maximum(i -> sum(abs, view(Pt, i, :) .- stat) / 2, axes(Pt, 1)) ≤ ε && return t
        Pt = Pt * m.P
    end
    return nothing
end
