function SoEwald2D_Es_pair(q_1::T, q_2::T, α::T, r_sq::T) where{T}
    return q_1 * q_2 * erfc(α * sqrt(r_sq)) / sqrt(r_sq) / T(2)
end

function SoEwald2D_Es_self(q::T, α::T) where{T}
    return - q^2 * α / sqrt(π)
end

function SoEwald2D_Es(interaction::SoEwald2DShortInteraction{T}, neighbor::CellList3D{T}, sys::MDSys{T}, info::SimulationInfo{T}) where{T}
    neighbor_list = neighbor.neighbor_list

    energy_short = zero(T)
    atoms = sys.atoms

    for (i, j, ρ) in neighbor_list
        id_i = info.particle_info[i].id
        id_j = info.particle_info[j].id
        coord_1, coord_2, r_sq = position_check3D(info.particle_info[i].position, info.particle_info[j].position, sys.boundary, interaction.r_c)
        if iszero(r_sq)
            nothing
        else
            q_1 = atoms[id_i].charge
            q_2 = atoms[id_j].charge
            energy_short += SoEwald2D_Es_pair(q_1, q_2, interaction.α, r_sq)
        end
    end

    for p_info in info.particle_info
        q = atoms[p_info.id].charge
        energy_short += SoEwald2D_Es_self(q, interaction.α)
    end

    return energy_short / (4π * interaction.ϵ_0)
end