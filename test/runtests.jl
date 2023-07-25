using SoEwald2D
using Test

@testset "SoEwald2D.jl" begin
    # Write your tests here.
    include("test_sort.jl")
    include("test_diff.jl")
    include("test_energy.jl")
end
