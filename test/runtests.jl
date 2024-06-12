using SpecialFunctions, ExTinyMD, EwaldSummations
using SoEwald2D
using Test

@testset "SoEwald2D.jl" begin
    include("soerfc.jl")
    include("energy.jl")
    include("force.jl")
    include("simulate.jl")
end
