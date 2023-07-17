module SoEwald2D

using SpecialFunctions, LinearAlgebra, Plots, Enzyme, GaussQuadrature, DelimitedFiles

include("types.jl")
include("soerfc.jl")
include("energy_long.jl")
include("force_long.jl")
include("direct_sum.jl")

end
