using SpecialFunctions, ExTinyMD, QuasiEwald, Distributed, StaticArrays, Random
using Test
@everywhere using SoEwald2D

@testset "SoEwald2D.jl" begin
    include("soerfc.jl")
    include("energy.jl")
    include("force.jl")
    include("simulate.jl")
    include("plan.jl")
    include("adapter.jl")
end
