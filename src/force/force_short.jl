# `coord_1`/`coord_2` are left untyped on purpose: `SVector{3,T}`,
# `NTuple{3,T}` and ExTinyMD's `Point{3,T}` all work with no conversion layer,
# since only `-`, `.-` and indexing are used. `sum(abs2, c1 .- c2)` replaces
# ExTinyMD's `dist2`, which is defined as exactly that for a `Point` pair.
#
# Sign: `ForwardDiff.derivative(energy, Δr)` is `dU/dr`, so `force` is
# `-dU/dr` and the returned vector is the FORCE on particle 1 -- repulsive
# (along `coord_1 - coord_2`) for like charges. That was already the
# convention here, unlike the long-range path; see `force!` below.
function SoEwald2D_Fs_pair(q_1::T, q_2::T, α::T, coord_1, coord_2) where{T}
    energy = r -> q_1 * q_2 * erfc(α * r) / r
    Δr = sqrt(sum(abs2, coord_1 .- coord_2))
    force = - ForwardDiff.derivative(energy, Δr)
    return force * (coord_1 - coord_2) / Δr
end

"""
    _short_pair_force(plan, poses, charges, i, j, r_c_sq) -> SVector{3,T}

Candidate-pair real-space force on particle `i` (the force on `j` is its
negation), zero unless the pair is inside the cutoff after the true
minimum-image correction. See [`_short_pair_energy`](@ref) for why the
`iszero(r_sq)` half of the guard is there and why zero rather than an error is
the right answer for a coincident pair -- here it is what stops a `0/0` NaN,
since the force is built from a unit vector.
"""
@inline function _short_pair_force(plan::SoEwald2DShortPlan{T}, poses, charges, i, j, r_c_sq::T) where{T}
    coord_1, coord_2, r_sq = _min_image_slab(poses[i], poses[j], plan.L)
    if r_sq ≥ r_c_sq || iszero(r_sq)
        return SVector{3, T}(zero(T), zero(T), zero(T))
    end
    return SoEwald2D_Fs_pair(charges[i], charges[j], plan.α, coord_1, coord_2)
end

"""
    SoEwald2D.force!(F, plan::SoEwald2DShortPlan, poses, charges; neighbor_list = nothing) -> F
    SoEwald2D.force(plan::SoEwald2DShortPlan, poses, charges; neighbor_list = nothing) -> Vector{SVector{3,T}}

Short-range (real-space) **force** from plain array-of-structs positions and
charges, written into `F` (filled, not accumulated into). Neither `poses` nor
`charges` is mutated. See [`SoEwald2D.energy`](@ref) for the `neighbor_list`
contract (candidate pairs only, separation always recomputed here) and for why
no `neighbor_list` means an `O(n_atoms^2)` pair loop.

Returns a force, i.e. `-∇U`, **not** the energy gradient, and no mass division
is applied -- that is the caller's job. The short-range path already used this
sign before the decoupling; the long-range one did not, see
[`SoEwald2D.force!`](@ref) for `SoEwald2DLongPlan`.

Not exported; call as `SoEwald2D.force!(...)`.
"""
function force!(F, plan::SoEwald2DShortPlan{T}, poses, charges; neighbor_list = nothing) where{T}
    n_atoms = plan.n_atoms
    r_c_sq = plan.r_c^2
    fill!(F, SVector{3, T}(zero(T), zero(T), zero(T)))

    if neighbor_list === nothing
        for i in 1:n_atoms, j in (i + 1):n_atoms
            F_ij = _short_pair_force(plan, poses, charges, i, j, r_c_sq) / (4π * plan.ϵ_0)
            F[i] += F_ij
            F[j] -= F_ij
        end
    else
        for pair in neighbor_list
            i, j = pair[1], pair[2]
            F_ij = _short_pair_force(plan, poses, charges, i, j, r_c_sq) / (4π * plan.ϵ_0)
            F[i] += F_ij
            F[j] -= F_ij
        end
    end

    return F
end

"Allocating form of [`SoEwald2D.force!`](@ref)."
function force(plan::SoEwald2DShortPlan{T}, poses, charges; kwargs...) where{T}
    F = [SVector{3, T}(zero(T), zero(T), zero(T)) for _ in 1:plan.n_atoms]
    return force!(F, plan, poses, charges; kwargs...)
end
