function force_sum_k0!(q::Vector{T}, z::Vector{T}, para::SoEwald2DLongInteraction{T}, soepara::SoePara{ComplexF64}, iterpara::IterPara, acceleration::Vector{Point{3, T}}) where{T}


    return nothing
end

function SoEwald2D_Fl!(interaction::SoEwald2DLongInteraction{T}, sys::MDSys, info::SimulationInfo{T}) where{T<:Number}

    revise_interaction!(interaction, sys, info)

    iterpara = interaction.iterpara
    soepara = interaction.soepara
    q = interaction.q
    x = interaction.x
    y = interaction.y
    z = interaction.z

    update_iterpara_z!(iterpara, z)

    force_sum_k0!(q, z, interaction, soepara, iterpara, info.acceleration)

    # info.acceleration .-= force_sum_k0(q, z, interaction, soepara, iterpara) / (4 * interaction.L[1] * interaction.L[2])

    for i in 1:size(interaction.k_set, 1)
        K = interaction.k_set[i]
        force_sum_k(K, q, x, y, z, interaction, soepara, iterpara, info.acceleration)
        # info.acceleration .+= force_sum_k(K, q, x, y, z, interaction, soepara, iterpara) / (8 * interaction.L[1] * interaction.L[2])
    end

    return nothing 
end