function SoEwald2D_Fs_pair(q_1::T, q_2::T, α::T, coord_1::Point{3, T}, coord_2::Point{3, T}) where{T}
    energy = r -> q_1 * q_2 * erfc(α * r) / r / T(2)
    Δr = sqrt(dist2(coord_1, coord_2))
    force = ForwardDiff.derivative(energy, Δr)
    return force * (coord_1 - coord_2) / Δr
end

function SoEwald2D_Fs!(interaction::SoEwald2DShortInteraction{T}, neighbor::CellList3D{T}, sys::MDSys{T}, info::SimulationInfo{T}) where{T}
    neighbor_list = neighbor.neighbor_list
    atoms = sys.atoms
    acceleration = info.acceleration

    for (i, j, ρ) in neighbor_list
        id_i = info.particle_info[i].id
        id_j = info.particle_info[j].id
        coord_1, coord_2, r_sq = position_check3D(info.particle_info[i].position, info.particle_info[j].position, sys.boundary, interaction.r_c)
        if iszero(r_sq)
            nothing
        else
            q_1 = atoms[id_i].charge
            q_2 = atoms[id_j].charge
            F_ij = SoEwald2D_Fs_pair(q_1, q_2, interaction.α, coord_1, coord_2)
            acceleration[id_i] += F_ij / (4π * interaction.ϵ_0)
            acceleration[id_j] -= F_ij / (4π * interaction.ϵ_0)
        end
    end

    return nothing
end