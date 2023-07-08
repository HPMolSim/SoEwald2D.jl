module SoEwald2D

using SpecialFunctions, LinearAlgebra, Plots, Enzyme, GaussQuadrature, DelimitedFiles

# Write your package code here.
include("types.jl")
include("soerfc.jl")
include("energy_long.jl")
include("force_long.jl")

end
