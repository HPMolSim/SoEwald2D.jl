using SoEwald2D
using Test, SpecialFunctions, ExTinyMD

@testset "SoEwald2D.jl" begin
    include("soerfc.jl")
    include("sorting.jl")
    include("energy.jl")
end
