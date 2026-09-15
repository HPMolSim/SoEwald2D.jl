module SoEwald2D

using SpecialFunctions, LinearAlgebra, Enzyme, GaussQuadrature, ForwardDiff, Distributions, Random, StatsBase, Distributed, StaticArrays

# `import ExTinyMD`, NOT `using ExTinyMD`. A blanket `using` brings ExTinyMD's
# exported `energy` into scope, and Julia refuses to let a module define its
# own `function energy(...)` while a `using`-imported binding of that name is
# visible -- even for a completely disjoint signature. Since this package now
# owns `SoEwald2D.energy`/`force`/`force!` (§4.2 of the decoupling design), the
# import has to be qualified. It disappears entirely in the next step, when
# src/adapter.jl moves to ext/SoEwald2DExTinyMDExt.jl.
import ExTinyMD

export SoePara, SoePara4, SoePara8, SoePara16, soerfc, soerf, soexp, soexp_mul_erfc
export AdPara, IterPara, revise_adpara!, update_iterpara_z!
export SoEwald2DShortPlan, SoEwald2DLongPlan
export SoEwald2DLongInteraction, SoEwald2DShortInteraction
export direct_sum, soe_direct_sum, diff_direct_sum

# `energy`, `force` and `force!` are declared here because their methods are
# spread across four included files, and they are deliberately NOT exported
# (§4.2): five packages in this family define an `energy`, so exporting them
# all would make the name ambiguous and error on first use under
# `using EwaldSummations, SoEwald2D`. Call them qualified:
# `SoEwald2D.energy(plan, poses, charges)`.
function energy end
function force end
function force! end

include("types.jl")

include("tools/soerfc.jl")
include("tools/direct_sum.jl")
include("tools/diff_direct_sum.jl")
include("tools/K_set.jl")

include("energy/energy_short.jl")
include("energy/energy_long.jl")

include("force/force_short.jl")
include("force/force_long.jl")

include("adapter.jl")

end
