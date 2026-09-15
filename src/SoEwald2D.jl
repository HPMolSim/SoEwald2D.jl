module SoEwald2D

using SpecialFunctions, LinearAlgebra, Enzyme, GaussQuadrature, ForwardDiff, Distributions, Random, StatsBase, Distributed, StaticArrays

# ExTinyMD is deliberately absent from this module. It is a [weakdeps] entry
# only; the ExTinyMD.AbstractInteraction wrappers and the ExTinyMD.energy /
# ExTinyMD.update_acceleration! methods all live in
# ext/SoEwald2DExTinyMDExt.jl. Everything below is plain arrays and numbers.
#
# (It was `import ExTinyMD`, never `using ExTinyMD`, for the one commit the
# adapter spent in src/: a blanket `using` brings ExTinyMD's exported `energy`
# into scope, and Julia refuses to let a module define its own
# `function energy(...)` while a `using`-imported binding of that name is
# visible, even for a disjoint signature.)

export SoePara, SoePara4, SoePara8, SoePara16, soerfc, soerf, soexp, soexp_mul_erfc
export AdPara, IterPara, revise_adpara!, update_iterpara_z!
export SoEwald2DShortPlan, SoEwald2DLongPlan
# SoEwald2DShortInteraction / SoEwald2DLongInteraction are
# ExTinyMD.AbstractInteraction wrappers that can only be *defined* once
# ExTinyMD exists: a struct's supertype is fixed where the struct is defined,
# and a struct definition cannot be dot-qualified at all -- `struct
# SoEwald2D.Foo <: ExTinyMD.AbstractInteraction ... end` is not legal Julia --
# so an extension cannot define a type INTO its parent package's namespace the
# way it defines a method into a parent's generic function. The wrapper types
# therefore live in the extension module's own namespace, and these two names
# are declared here as plain dispatcher functions -- exported, so
# `using SoEwald2D, ExTinyMD; SoEwald2DShortInteraction(...)` keeps working
# exactly as before -- that hand off to the extension's real constructors via
# `Base.get_extension`.
#
# The cost, and it was checked before choosing this over renaming: the
# preserved name is a FUNCTION, not a type. Construction works unchanged; every
# type-position use breaks. Grepped this repository and FastSpecSoG.jl (the only
# other package in the family that could plausibly depend on these names):
# FastSpecSoG mentions SoEwald2D nowhere at all, and inside this repository
# every type-position use was an internal `::SoEwald2DLongInteraction{T}`
# signature in src/, all of which now take the plan type instead. Outside src/,
# the tests and the README only ever CONSTRUCT. So the pattern is free here.
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

"""
    SoEwald2DShortInteraction(ϵ_0, L, s, α, n_atoms, r_c)
    SoEwald2DLongInteraction(ϵ_0, L, s, α, n_atoms, k_c, soepara; rbm = false, rbm_p = 0, parallel = true, rng = MersenneTwister(123))

`ExTinyMD.AbstractInteraction` wrappers around a [`SoEwald2DShortPlan`](@ref) /
[`SoEwald2DLongPlan`](@ref), for use in `sys.interactions`. Defined by
`ext/SoEwald2DExTinyMDExt.jl`, loaded automatically once `using ExTinyMD` has
also been done -- calling these before that raises an informative error rather
than a cryptic `UndefVarError`. For standalone use (no ExTinyMD, no MDSys),
construct a plan directly and query it with
`SoEwald2D.energy`/`force`/`force!`.

!!! warning "Breaking change in 0.2.0: these two names are functions, not types"
    Before this package was decoupled from ExTinyMD, each was a `struct`. They
    are now *dispatcher functions* that forward to the real constructors in the
    extension, so:

    * **Construction works unchanged.** `SoEwald2DShortInteraction(ϵ_0, L, s,
      α, n_atoms, r_c)` returns exactly the wrapper it always did, and it still
      `isa ExTinyMD.AbstractInteraction`.
    * **Type-position uses do not work.** `x isa SoEwald2DShortInteraction`, an
      `::SoEwald2DLongInteraction` annotation, a
      `Vector{SoEwald2DShortInteraction}` element type, and dispatching a
      method on one all now raise a `TypeError`.

    If you need the type itself, reach into the extension module:

    ```julia
    using SoEwald2D, ExTinyMD
    ext = Base.get_extension(SoEwald2D, :SoEwald2DExTinyMDExt)
    x isa ext.SoEwald2DShortInteraction        # works
    ```

    The type cannot be re-exported from this module under the same name,
    because the exported name is what makes the constructor call resolve
    without ExTinyMD being a hard dependency.
"""
SoEwald2DShortInteraction, SoEwald2DLongInteraction

# ExTinyMD's PkgId, for telling "ExTinyMD was never loaded" apart from
# "ExTinyMD is loaded but the extension did not come up". `Base.get_extension`
# returns `nothing` in both cases, and reporting the second as the first sends
# the user off to run a `using ExTinyMD` they have already run, while the
# actual precompilation error has scrolled off the screen. The second case is
# the one you actually hit while developing the extension.
const _EXTINYMD_PKGID = Base.PkgId(Base.UUID("fec76197-d59f-46dd-a0ed-76a83c21f7aa"), "ExTinyMD")

for name in (:SoEwald2DShortInteraction, :SoEwald2DLongInteraction)
    @eval function $name(args...; kwargs...)
        ext = Base.get_extension(SoEwald2D, :SoEwald2DExTinyMDExt)
        if ext === nothing
            if haskey(Base.loaded_modules, _EXTINYMD_PKGID)
                error($(string(name)) * ": ExTinyMD is loaded, but SoEwald2D's " *
                      "SoEwald2DExTinyMDExt extension failed to load -- so `using ExTinyMD` " *
                      "is not what is missing. This is almost always a precompilation error " *
                      "inside the extension, reported as an `Error: Error during loading of " *
                      "extension SoEwald2DExTinyMDExt of SoEwald2D` warning that has since " *
                      "scrolled past. Run `Base.retry_load_extensions()` to reproduce it, and " *
                      "check for a version conflict between the loaded ExTinyMD and this " *
                      "package's [compat] bound.")
            else
                error($(string(name)) * " requires ExTinyMD to be loaded (`using ExTinyMD`) -- " *
                      "it is defined by the SoEwald2DExTinyMDExt package extension. For standalone " *
                      "use, construct a SoEwald2DShortPlan/SoEwald2DLongPlan directly instead and " *
                      "query it with SoEwald2D.energy / SoEwald2D.force / SoEwald2D.force!.")
            end
        end
        return getfield(ext, $(QuoteNode(name)))(args...; kwargs...)
    end
end

include("types.jl")

include("tools/soerfc.jl")
include("tools/direct_sum.jl")
include("tools/diff_direct_sum.jl")
include("tools/K_set.jl")

include("energy/energy_short.jl")
include("energy/energy_long.jl")

include("force/force_short.jl")
include("force/force_long.jl")

end
