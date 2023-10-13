function ExTinyMD.update_acceleration!(interaction::SoEwald2DShortInteraction{T}, neighborfinder::T_NEIGHBOR, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, T_NEIGHBOR<:ExTinyMD.AbstractNeighborFinder}
    update_finder!(neighborfinder, info)
    SoEwald2D_Fs!(interaction, neighborfinder, sys, info)
    return nothing
end

function ExTinyMD.update_acceleration!(interaction::SoEwald2DLongInteraction{T}, neighborfinder::T_NEIGHBOR, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, T_NEIGHBOR<:ExTinyMD.AbstractNeighborFinder}

    SoEwald2D_Fl!(interaction, sys, info)
    return nothing
end