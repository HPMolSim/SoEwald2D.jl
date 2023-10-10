module SoEwald2D

using SpecialFunctions, LinearAlgebra, Enzyme, GaussQuadrature, ExTinyMD, ForwardDiff

export SoePara, SoEwald2DPara, AdPara, IterPara, revise_adpara!, update_iterpara!, revise_interaction!, update_iterpara_z!, update_iterpara_m!
export SoEwald2DLongInteraction, SoEwald2DShortInteraction, SoEwald2D_Fs!, SoEwald2D_Fl!, SoEwald2D_El, SoEwald2D_Es

include("types.jl")

include("tools/soerfc.jl")
include("tools/direct_sum.jl")

include("energy/energy.jl")
include("energy/energy_short.jl")
include("energy/energy_long.jl")

include("force/force.jl")
include("force/force_short.jl")
include("force/force_long.jl")

end
