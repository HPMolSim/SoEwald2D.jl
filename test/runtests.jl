using SoEwald2D
using SpecialFunctions, ExTinyMD, QuasiEwald
using Test

@testset "SoEwald2D.jl" begin
    include("soerfc.jl")
    include("sorting.jl")
    include("energy.jl")
end
