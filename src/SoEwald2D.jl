module SoEwald2D

using SpecialFunctions, LinearAlgebra, Enzyme, GaussQuadrature

include("types.jl")
include("soerfc.jl")
include("energy_long.jl")
include("energy_long_split.jl")
include("force_long.jl")
include("direct_sum.jl")

end
