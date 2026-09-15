# `coord_1`/`coord_2` are left untyped on purpose: SVector{3,T}, NTuple{3,T} and
# ExTinyMD's Point{3,T} all work with no conversion layer, since only `-` and
# indexing are used. `sum(abs2, c1 .- c2)` replaces ExTinyMD's `dist2`, which is
# defined as exactly that for a Point pair.
function SoEwald2D_Fs_pair(q_1::T, q_2::T, α::T, coord_1, coord_2) where{T}
    energy = r -> q_1 * q_2 * erfc(α * r) / r
    Δr = sqrt(sum(abs2, coord_1 .- coord_2))
    force = - ForwardDiff.derivative(energy, Δr)
    return force * (coord_1 - coord_2) / Δr
end

function SoEwald2D_Fs!(interaction::SoEwald2DShortInteraction{T}, neighbor::CellList3D{T}, sys::MDSys{T}, info::SimulationInfo{T}) where{T}
    neighbor_list = neighbor.neighbor_list
    atoms = sys.atoms

    for (i, j, ρ) in neighbor_list
        id_i = info.particle_info[i].id
        id_j = info.particle_info[j].id
        coord_1, coord_2, r_sq = _min_image_slab(info.particle_info[i].position, info.particle_info[j].position, interaction.L)
        if r_sq ≥ interaction.r_c^2 || iszero(r_sq)
            nothing
        else
            q_1 = atoms[id_i].charge
            q_2 = atoms[id_j].charge
            F_ij = SoEwald2D_Fs_pair(q_1, q_2, interaction.α, coord_1, coord_2)
            F_ij_scaled = F_ij / (4π * interaction.ϵ_0)
            info.particle_info[i].acceleration += Point(F_ij_scaled[1], F_ij_scaled[2], F_ij_scaled[3])
            info.particle_info[j].acceleration -= Point(F_ij_scaled[1], F_ij_scaled[2], F_ij_scaled[3])
        end
    end

    return nothing
end