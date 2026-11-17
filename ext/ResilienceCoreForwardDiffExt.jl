module ResilienceCoreForwardDiffExt

# Automatic differentiation for linearise: the Jacobian of a SciML right-hand side.

using ForwardDiff: ForwardDiff
using ResilienceCore: ResilienceCore

function ResilienceCore.autodiff_jacobian(f, u::AbstractVector{<:Real}, inplace::Bool)
    x = Vector{Float64}(u)
    if inplace
        return Matrix{Float64}(ForwardDiff.jacobian(f, similar(x), x))
    end
    return Matrix{Float64}(ForwardDiff.jacobian(f, x))
end

end
