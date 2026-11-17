module ResilienceCoreSciMLBaseExt

# linearise for SciML problems: the Jacobian at an equilibrium (ODE, steady state) or a fixed
# point (discrete map). Documented on the `linearise` binding in ResilienceCore.

using LinearAlgebra: norm
using ResilienceCore: ResilienceCore, ContinuousJacobian, DiscreteJacobian, autodiff_jacobian
using SciMLBase: SciMLBase, ODEProblem, SteadyStateProblem, DiscreteProblem, isinplace

# f(u) of the problem, in place or not
function rhs(prob, u, p, t)
    if isinplace(prob)
        du = similar(u, float(eltype(u)))
        prob.f(du, u, p, t)
        return du
    end
    return prob.f(u, p, t)
end

analytic_jac(f) = hasproperty(f, :jac) ? getproperty(f, :jac) : nothing

function problem_jacobian(prob, u, p, t)
    f = prob.f
    jac = analytic_jac(f)
    n = length(u)
    if jac !== nothing
        if isinplace(prob)
            J = zeros(n, n)
            jac(J, u, p, t)
            return J
        end
        return Matrix{Float64}(jac(u, p, t))
    end
    if isinplace(prob)
        return autodiff_jacobian((du, x) -> f(du, x, p, t), u, true)
    end
    return autodiff_jacobian(x -> f(x, p, t), u, false)
end

function check_point(residual, u, tol, what)
    r = norm(residual, Inf)
    bound = tol * max(1.0, norm(u, Inf))
    r <= bound ||
        throw(ArgumentError("u is not $what of the problem: the residual is $r, above the tolerance $bound. " *
                            "Find the point first (e.g. with a SteadyStateProblem), raise `tol`, or pass " *
                            "`check = false`."))
    return nothing
end

function continuous(prob, u, p, t, tol, check, equilibrium)
    u isa AbstractVector || throw(ArgumentError("the state must be a vector, got $(typeof(u))"))
    check && check_point(rhs(prob, u, p, t), u, tol, "an equilibrium (f(u) = 0)")
    return ContinuousJacobian(problem_jacobian(prob, u, p, t), equilibrium)
end

function ResilienceCore.linearise(prob::ODEProblem, u = prob.u0; p = prob.p, t = first(prob.tspan),
                                  tol = 1e-8, check = true, equilibrium = u)
    return continuous(prob, u, p, t, tol, check, equilibrium)
end

function ResilienceCore.linearise(prob::SteadyStateProblem, u = prob.u0; p = prob.p, t = 0.0, tol = 1e-8,
                                  check = true, equilibrium = u)
    return continuous(prob, u, p, t, tol, check, equilibrium)
end

function ResilienceCore.linearise(prob::DiscreteProblem, u = prob.u0; p = prob.p, t = first(prob.tspan),
                                  tol = 1e-8, check = true, equilibrium = u)
    u isa AbstractVector || throw(ArgumentError("the state must be a vector, got $(typeof(u))"))
    check && check_point(rhs(prob, u, p, t) .- u, u, tol, "a fixed point (f(u) = u)")
    return DiscreteJacobian(problem_jacobian(prob, u, p, t), equilibrium)
end

end
