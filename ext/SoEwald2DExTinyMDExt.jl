module SoEwald2DExTinyMDExt

# Bridge between ExTinyMD's MD loop and SoEwald2D's framework-free plans.
#
# ## Why the interaction types live here, not in src/
#
# `MDSys`'s constructor requires `interactions::Vector{T_INTERACTION}` with
# `T_INTERACTION <: Tuple{ExTinyMD.AbstractInteraction, ExTinyMD.AbstractNeighborFinder}`.
# A struct's supertype is fixed where the struct is defined -- Julia has no
# mechanism for an extension, loaded later and conditionally, to retroactively
# add a supertype to an already-compiled type. `SoEwald2D`'s src/ does not
# depend on ExTinyMD at all (a weak dependency only), so nothing defined there
# can ever be a subtype of `ExTinyMD.AbstractInteraction` -- not "unless you
# remember an annotation", but structurally, in every build of the package.
#
# `SoEwald2DShortPlan`/`SoEwald2DLongPlan` (src/types.jl) are the framework-free
# core: parameters plus, for the long plan, its structure-of-arrays solver
# scratch, and no notion of mass at all. The two structs below are thin
# `ExTinyMD.AbstractInteraction` wrappers around one of those, holding only the
# MD-side scratch (gathered positions/charges, a force buffer) that a plan has
# no business owning. `SoEwald2D.jl`'s main module declares
# `SoEwald2DShortInteraction`/`SoEwald2DLongInteraction` as dispatcher functions
# that call through to the constructors below via `Base.get_extension` once this
# extension has loaded, so `using SoEwald2D, ExTinyMD;
# SoEwald2DShortInteraction(...)` keeps working exactly as it did before this
# package was decoupled.

using SoEwald2D, ExTinyMD, StaticArrays

# ----------------------------------------------------------------------------
# Gather helpers. Index convention (matching ExTinyMD's own adapter,
# ../ExTinyMD.jl/src/interactions/electrostatics/adapter.jl):
# `info.particle_info` is indexed by storage SLOT, `sys.atoms` by particle ID.
# Positions are read in slot order and charges/masses gathered to match, so
# plan-layer index `i` consistently means "slot i" and forces come back in slot
# order. In stock ExTinyMD `particle_info[i].id == i`, so a gather that confuses
# the two looks correct forever -- test/adapter.jl permutes the mapping and
# gives every id a distinct mass and charge so that it cannot.
# ----------------------------------------------------------------------------

"Gather positions in slot order into `buf`, as SVector{3,T} (the plan's canonical AoS element)."
function gather_positions!(buf::Vector{SVector{3, T}}, info::ExTinyMD.SimulationInfo{T}) where {T}
    @inbounds for i in eachindex(info.particle_info)
        p = info.particle_info[i].position
        buf[i] = SVector{3, T}(p[1], p[2], p[3])
    end
    return buf
end

"Gather charges in slot order into `buf`, honouring the id/slot indirection."
function gather_charges!(buf::Vector{T}, sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where {T}
    @inbounds for i in eachindex(info.particle_info)
        buf[i] = sys.atoms[info.particle_info[i].id].charge
    end
    return buf
end

# Which finders the short-range wrapper accepts, and why anything else is a
# hard error rather than a silent `f.neighbor_list`.
#
# SoEwald2D's real-space sum is `erfc(α·r)/r` in the FULL three-dimensional
# separation, so its candidate pairs have to be selected by 3-D distance.
# `CellList3D` and `CellListDir3D` do exactly that. `NoNeighborFinder` carries
# no list at all, so the plan falls back to its own O(n^2) pair loop, which is
# always correct and merely slower.
#
# Everything else must be a hard error, NOT a silent `f.neighbor_list`. A
# quasi-2D finder (`CellListQ2D`, `CellListDirQ2D`) has a `neighbor_list` field
# of the right shape, so an untyped fallback accepts it happily -- and then
# feeds this sum a candidate set chosen by in-plane distance, which is a
# different set: it omits nothing the 3-D cutoff needs but includes tall
# near-columnar pairs the plan then has to reject one by one, so a finder
# mismatch shows up as a silent performance cliff rather than an error, and a
# *smaller* Q2D cutoff would silently drop pairs outright. Before the
# decoupling, `SoEwald2D_Es`/`SoEwald2D_Fs!` were annotated `::CellList3D{T}`
# and anything else was a `MethodError`; this restores that loud failure with a
# message that says what to use instead.
_finder_list(::ExTinyMD.NoNeighborFinder) = nothing
_finder_list(f::Union{ExTinyMD.CellList3D, ExTinyMD.CellListDir3D}) = f.neighbor_list
_finder_list(f) = throw(ArgumentError(
    "SoEwald2D's short-range interaction needs a 3-D neighbour finder, but got " *
    "a $(typeof(f)). Use ExTinyMD.CellList3D or ExTinyMD.CellListDir3D (pairs " *
    "selected by full three-dimensional distance, which is what the real-space " *
    "kernel erfc(α·r)/r needs), or ExTinyMD.NoNeighborFinder to fall back to " *
    "the plan's own O(n^2) pair loop. A quasi-2D finder (CellListQ2D, " *
    "CellListDirQ2D) is refused deliberately and not merely unimplemented: its " *
    "candidate list is selected by in-plane distance, so it does not match this " *
    "sum's cutoff and a mismatch would be silent."))

# `+=`, never `=`: `update_acceleration!` is called once per interaction per
# step and each one must add to what the others already put there.
#
# The mass division belongs here, not in the plan: a framework-free solver
# returns a force. Bitwise-equivalent to what the long-range path did before
# (`acceleration -= Point(Fx, Fy, Fz) / mass`), because ExTinyMD's `Point` `/`
# and `-` are plain componentwise operations and IEEE-754 makes
# `a - (f/m) == a + ((-f)/m)` exactly -- negation is exact and division is
# correctly rounded. Verified against a captured baseline, not just argued.
#
# It is NOT equivalent for the short-range path, and that is a deliberate fix:
# the pre-decoupling `SoEwald2D_Fs!` did `acceleration += F_ij / (4π ϵ_0)` with
# no mass division at all, so it added a force to an acceleration. Every test
# in the old suite used mass = 1.0, which is why it went unnoticed. Reported,
# not absorbed.
function _accumulate_acceleration!(F, sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where {T}
    @inbounds for i in eachindex(info.particle_info)
        m = sys.atoms[info.particle_info[i].id].mass
        f = F[i]
        info.particle_info[i].acceleration += ExTinyMD.Point(f[1] / m, f[2] / m, f[3] / m)
    end
    return nothing
end

# ----------------------------------------------------------------------------
# Short-range wrapper
# ----------------------------------------------------------------------------

"""
    SoEwald2DShortInteraction(ϵ_0, L, s, α, n_atoms, r_c)
    SoEwald2DShortInteraction(plan::SoEwald2D.SoEwald2DShortPlan)

`ExTinyMD.AbstractInteraction` wrapper around a
[`SoEwald2D.SoEwald2DShortPlan`](@ref). Construct exactly as the
pre-decoupling `SoEwald2DShortInteraction` was constructed; place
`(interaction, finder)` in `sys.interactions` as before, where `finder` is an
`ExTinyMD.CellList3D`, `CellListDir3D` or `NoNeighborFinder`.
"""
struct SoEwald2DShortInteraction{P, T} <: ExTinyMD.AbstractInteraction
    plan::P
    pos_scratch::Vector{SVector{3, T}}
    charge_scratch::Vector{T}
    force_buffer::Vector{SVector{3, T}}
end

function SoEwald2DShortInteraction(plan::SoEwald2D.SoEwald2DShortPlan{T}) where {T}
    n = plan.n_atoms
    return SoEwald2DShortInteraction{typeof(plan), T}(
        plan, Vector{SVector{3, T}}(undef, n), Vector{T}(undef, n), Vector{SVector{3, T}}(undef, n))
end

SoEwald2DShortInteraction(ϵ_0, L, s, α, n_atoms, r_c) =
    SoEwald2DShortInteraction(SoEwald2D.SoEwald2DShortPlan(ϵ_0, L, s, α, n_atoms, r_c))

function ExTinyMD.energy(inter::SoEwald2DShortInteraction, neighborfinder,
                         sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where {T}
    ExTinyMD.update_finder!(neighborfinder, info)
    poses = gather_positions!(inter.pos_scratch, info)
    charges = gather_charges!(inter.charge_scratch, sys, info)
    return SoEwald2D.energy(inter.plan, poses, charges; neighbor_list = _finder_list(neighborfinder))
end

function ExTinyMD.update_acceleration!(inter::SoEwald2DShortInteraction, neighborfinder,
                                       sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where {T}
    ExTinyMD.update_finder!(neighborfinder, info)
    poses = gather_positions!(inter.pos_scratch, info)
    charges = gather_charges!(inter.charge_scratch, sys, info)
    F = SoEwald2D.force!(inter.force_buffer, inter.plan, poses, charges;
                         neighbor_list = _finder_list(neighborfinder))
    _accumulate_acceleration!(F, sys, info)
    return nothing
end

# ----------------------------------------------------------------------------
# Long-range wrapper
# ----------------------------------------------------------------------------

"""
    SoEwald2DLongInteraction(ϵ_0, L, s, α, n_atoms, k_c, soepara; rbm = false, rbm_p = 0, parallel = true, rng = MersenneTwister(123))
    SoEwald2DLongInteraction(plan::SoEwald2D.SoEwald2DLongPlan)

`ExTinyMD.AbstractInteraction` wrapper around a
[`SoEwald2D.SoEwald2DLongPlan`](@ref). Construct exactly as the
pre-decoupling `SoEwald2DLongInteraction` was constructed; place
`(interaction, finder)` in `sys.interactions` as before.

The neighbour finder is accepted and ignored -- the reciprocal-space sum runs
over every particle and consumes no candidate list -- so `NoNeighborFinder()`
is the natural choice and nothing needs restricting here, unlike the
short-range wrapper.
"""
struct SoEwald2DLongInteraction{P, T} <: ExTinyMD.AbstractInteraction
    plan::P
    pos_scratch::Vector{SVector{3, T}}
    charge_scratch::Vector{T}
    force_buffer::Vector{SVector{3, T}}
end

function SoEwald2DLongInteraction(plan::SoEwald2D.SoEwald2DLongPlan{T}) where {T}
    n = plan.n_atoms
    return SoEwald2DLongInteraction{typeof(plan), T}(
        plan, Vector{SVector{3, T}}(undef, n), Vector{T}(undef, n), Vector{SVector{3, T}}(undef, n))
end

SoEwald2DLongInteraction(ϵ_0, L, s, α, n_atoms, k_c, soepara; kwargs...) =
    SoEwald2DLongInteraction(SoEwald2D.SoEwald2DLongPlan(ϵ_0, L, s, α, n_atoms, k_c, soepara; kwargs...))

function ExTinyMD.energy(inter::SoEwald2DLongInteraction, neighborfinder,
                         sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where {T}
    poses = gather_positions!(inter.pos_scratch, info)
    charges = gather_charges!(inter.charge_scratch, sys, info)
    return SoEwald2D.energy(inter.plan, poses, charges)
end

function ExTinyMD.update_acceleration!(inter::SoEwald2DLongInteraction, neighborfinder,
                                       sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where {T}
    poses = gather_positions!(inter.pos_scratch, info)
    charges = gather_charges!(inter.charge_scratch, sys, info)
    F = SoEwald2D.force!(inter.force_buffer, inter.plan, poses, charges)
    _accumulate_acceleration!(F, sys, info)
    return nothing
end

end
