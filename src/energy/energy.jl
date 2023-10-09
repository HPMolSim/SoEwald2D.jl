function ExTinyMD.energy(interaction::SoEwald2DShortInteraction{T}, neighborfinder::T_NEIGHBOR, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, T_NEIGHBOR<:ExTinyMD.AbstractNeighborFinder}
    update_finder!(neighborfinder, info)
    return QuasiEwald_Es(interaction, neighborfinder, sys, info)
end

function ExTinyMD.energy(interaction::SoEwald2DLongInteraction{T}, neighborfinder::T_NEIGHBOR, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, T_NEIGHBOR<:ExTinyMD.AbstractNeighborFinder}
    
    return QuasiEwald_El(interaction, sys, info)
end