using ResilienceCore
using Test
using LinearAlgebra
using StructuredPopulationCore: StructuredPopulationCore

@testset "ResilienceCore" begin
    include("test_dispatch.jl")
    include("test_consistency.jl")
    include("test_projection.jl")
    include("test_markov_periodic.jl")
    include("test_community.jl")
end
