module SoEwald2D

using SpecialFunctions, LinearAlgebra, Enzyme, GaussQuadrature, ExTinyMD, ForwardDiff, Distributions, Random, StatsBase, Distributed

export SoePara, SoePara4, SoePara8, SoePara16, soerfc, soerf, soexp, soexp_mul_erfc
export AdPara, IterPara, revise_adpara!, update_iterpara!, revise_interaction!, update_iterpara_z!
export SoEwald2D_init, SoEwald2DLongInteraction, SoEwald2DShortInteraction, SoEwald2D_Fs!, SoEwald2D_Fl!, SoEwald2D_El, SoEwald2D_Es
export direct_sum, soe_direct_sum, diff_direct_sum, ad_energy_sum!

include("types.jl")

include("tools/soerfc.jl")
include("tools/direct_sum.jl")
include("tools/diff_direct_sum.jl")
include("tools/K_set.jl")

include("energy/energy.jl")
include("energy/energy_short.jl")
include("energy/energy_long.jl")

include("force/force.jl")
include("force/force_short.jl")
include("force/force_long.jl")

end
