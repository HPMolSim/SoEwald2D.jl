using SpecialFunctions, ExTinyMD, QuasiEwald, Distributed
using Test
@everywhere using SoEwald2D

@testset "SoEwald2D.jl" begin
    include("soerfc.jl")
    include("energy.jl")
    include("force.jl")
    include("simulate.jl")
end
