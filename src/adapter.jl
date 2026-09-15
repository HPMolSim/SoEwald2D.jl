# ============================================================================
# The ExTinyMD bridge. Everything in this file moves verbatim into
# ext/SoEwald2DExTinyMDExt.jl in the next step -- it is the last
# ExTinyMD-coupled code in src/, and it is kept here for exactly one commit so
# that the core can be swapped out underneath it with the suite staying green.
#
# Index convention, matching ExTinyMD's own adapter
# (../ExTinyMD.jl/src/interactions/electrostatics/adapter.jl):
# `info.particle_info` is indexed by storage SLOT, `sys.atoms` by particle ID.
# Positions are read in slot order and charges/masses gathered to match, so
# plan-layer index `i` consistently means "slot i" and forces come back in slot
# order.
# ============================================================================

"Gather positions in slot order into `buf`, as SVector{3,T} (the plan's canonical AoS element)."
function gather_positions!(buf::Vector{SVector{3, T}}, info::ExTinyMD.SimulationInfo{T}) where{T}
    @inbounds for i in eachindex(info.particle_info)
        p = info.particle_info[i].position
        buf[i] = SVector{3, T}(p[1], p[2], p[3])
    end
    return buf
end

"Gather charges in slot order into `buf`, honouring the id/slot indirection."
function gather_charges!(buf::Vector{T}, sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where{T}
    @inbounds for i in eachindex(info.particle_info)
        buf[i] = sys.atoms[info.particle_info[i].id].charge
    end
    return buf
end

# Which finders the short-range wrapper accepts, and why anything else is a
# hard error rather than a silent `f.neighbor_list`.
#
# SoEwald2D's real-space sum needs candidate pairs selected by the FULL 3-D
# distance -- its kernel is `erfc(α·r)/r` in the true separation. `CellList3D`
# and `CellListDir3D` do exactly that; `NoNeighborFinder` carries no list at
# all, so the plan falls back to its own O(n^2) pair loop, which is always
# correct and merely slower. A quasi-2D finder (`CellListQ2D`,
# `CellListDirQ2D`) selects on the in-plane distance instead, so it would
# *over*-supply candidates -- harmless in itself, since the plan recomputes
# every separation -- but it is refused anyway so that a finder mismatch is
# loud rather than a silent performance cliff, and so that the error message
# can say which finder this interaction wants. Before the decoupling,
# `SoEwald2D_Es`/`SoEwald2D_Fs!` were annotated `::CellList3D{T}` and anything
# else was a `MethodError`; this keeps that loud failure with a better message.
_finder_list(::ExTinyMD.NoNeighborFinder) = nothing
_finder_list(f::Union{ExTinyMD.CellList3D, ExTinyMD.CellListDir3D}) = f.neighbor_list
_finder_list(f) = throw(ArgumentError(
    "SoEwald2D's short-range interaction needs a 3-D neighbour finder, but got " *
    "a $(typeof(f)). Use ExTinyMD.CellList3D or ExTinyMD.CellListDir3D (pairs " *
    "selected by full three-dimensional distance, which is what the real-space " *
    "kernel erfc(α·r)/r needs), or ExTinyMD.NoNeighborFinder to fall back to the " *
    "plan's own O(n^2) pair loop."))

function ExTinyMD.energy(inter::SoEwald2DShortInteraction{T}, neighborfinder, sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where{T}
    ExTinyMD.update_finder!(neighborfinder, info)
    poses = gather_positions!(inter.pos_scratch, info)
    charges = gather_charges!(inter.charge_scratch, sys, info)
    return energy(inter.plan, poses, charges; neighbor_list = _finder_list(neighborfinder))
end

function ExTinyMD.energy(inter::SoEwald2DLongInteraction{T}, neighborfinder, sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where{T}
    poses = gather_positions!(inter.pos_scratch, info)
    charges = gather_charges!(inter.charge_scratch, sys, info)
    return energy(inter.plan, poses, charges)
end

# `+=`, never `=`: `update_acceleration!` is called once per interaction per
# step and each one must add to what the others already put there. The mass
# division belongs here, not in the plan -- a framework-free solver returns a
# force.
function _accumulate_acceleration!(F, sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where{T}
    @inbounds for i in eachindex(info.particle_info)
        m = sys.atoms[info.particle_info[i].id].mass
        f = F[i]
        info.particle_info[i].acceleration += ExTinyMD.Point(f[1] / m, f[2] / m, f[3] / m)
    end
    return nothing
end

function ExTinyMD.update_acceleration!(inter::SoEwald2DShortInteraction{T}, neighborfinder, sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where{T}
    ExTinyMD.update_finder!(neighborfinder, info)
    poses = gather_positions!(inter.pos_scratch, info)
    charges = gather_charges!(inter.charge_scratch, sys, info)
    F = force!(inter.force_buffer, inter.plan, poses, charges; neighbor_list = _finder_list(neighborfinder))
    _accumulate_acceleration!(F, sys, info)
    return nothing
end

function ExTinyMD.update_acceleration!(inter::SoEwald2DLongInteraction{T}, neighborfinder, sys::ExTinyMD.MDSys{T}, info::ExTinyMD.SimulationInfo{T}) where{T}
    poses = gather_positions!(inter.pos_scratch, info)
    charges = gather_charges!(inter.charge_scratch, sys, info)
    F = force!(inter.force_buffer, inter.plan, poses, charges)
    _accumulate_acceleration!(F, sys, info)
    return nothing
end
