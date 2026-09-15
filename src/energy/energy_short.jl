function SoEwald2D_Es_pair(q_1::T, q_2::T, α::T, r_sq::T) where{T}
    return q_1 * q_2 * erfc(α * sqrt(r_sq)) / sqrt(r_sq)
end

function SoEwald2D_Es_self(q::T, α::T) where{T}
    return - q^2 * α / sqrt(π)
end

"""
    _short_pair_energy(plan, poses, charges, i, j, r_c_sq) -> T

Candidate-pair real-space energy, zero unless the pair is inside the cutoff
after the true minimum-image correction.

Two guards, not one. `r_sq ≥ r_c_sq` is the cutoff test that replaces
`position_check3D`'s all-zero sentinel triple. `iszero(r_sq)` is the *other*
job that sentinel was quietly doing: it also skipped coincident pairs. Keeping
it is deliberate --

  * `erfc(α·r)/r` diverges as `r → 0`, so no finite value is "the limit"; and
    the matching force is built from a unit vector, so dropping this guard puts
    a `0/0` NaN in the force array (the exact bug that shipped in QuasiEwald
    when this substitution was made there, where it fired for every lattice
    initialisation).
  * `r == 0` requires all three separations to vanish after the x/y wrap, i.e.
    two particles at the same point. That is a configuration error, not a
    physical state, and one such pair should not turn every other particle's
    force into a NaN.
  * Returning zero reproduces the pre-decoupling behaviour bit for bit, which
    keeps this refactor free of numerical change.

An `ArgumentError` was the alternative and was rejected: a coincident pair is
reachable from a perfectly ordinary lattice generator (two sites whose x/y
separation is an exact multiple of `Lx`/`Ly` at equal `z`), and turning that
into a hard failure would be a behaviour change on top of a refactor that is
otherwise numerically inert.
"""
@inline function _short_pair_energy(plan::SoEwald2DShortPlan{T}, poses, charges, i, j, r_c_sq::T) where{T}
    _, _, r_sq = _min_image_slab(poses[i], poses[j], plan.L)
    if r_sq ≥ r_c_sq || iszero(r_sq)
        return zero(T)
    end
    return SoEwald2D_Es_pair(charges[i], charges[j], plan.α, r_sq)
end

"""
    SoEwald2D.energy(plan::SoEwald2DShortPlan, poses, charges; neighbor_list = nothing) -> T

Short-range (real-space) energy from plain array-of-structs positions and
charges. No ExTinyMD type is constructed and neither argument is mutated.

Pass `neighbor_list` (an iterable of `(i, j, ...)` candidate pairs, e.g. the
`neighbor_list` field of an ExTinyMD `CellList3D`) to reuse a list maintained
elsewhere. Each entry is treated as a *candidate only*: the separation is
always recomputed here from `poses` via [`_min_image_slab`](@ref), and a
supplied list's own reported distance is ignored. That discipline is not
pedantry -- Phase 1 of this project produced a wrong-signed `Ewald2D` energy
(+0.0238 against a true −0.1539) precisely by trusting a neighbour list's `r`,
which was an in-plane distance from a 2-D finder.

With no `neighbor_list` every pair is tested directly, `O(n_atoms^2)`; this
plan owns no cell list of its own. With one, the result depends on the list's
order only through floating-point summation order.

`energy` is deliberately **not exported** -- call it as
`SoEwald2D.energy(...)`. Five packages in this family define an `energy`, and
exporting them all would make the name ambiguous on `using`.
"""
function energy(plan::SoEwald2DShortPlan{T}, poses, charges; neighbor_list = nothing) where{T}
    n_atoms = plan.n_atoms
    r_c_sq = plan.r_c^2
    energy_short = zero(T)

    if neighbor_list === nothing
        for i in 1:n_atoms, j in (i + 1):n_atoms
            energy_short += _short_pair_energy(plan, poses, charges, i, j, r_c_sq)
        end
    else
        for pair in neighbor_list
            energy_short += _short_pair_energy(plan, poses, charges, pair[1], pair[2], r_c_sq)
        end
    end

    for i in 1:n_atoms
        energy_short += SoEwald2D_Es_self(charges[i], plan.α)
    end

    return energy_short / (4π * plan.ϵ_0)
end
