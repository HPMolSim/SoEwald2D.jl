Base.real(x::Point) = Point(real.(x.coo))

function ad_energy_sum!(interaction::SoEwald2DLongInteraction{T}) where{T}

    iterpara = interaction.iterpara
    soepara = interaction.soepara
    adpara = interaction.adpara
    q = interaction.q
    x = interaction.x
    y = interaction.y
    z = interaction.z
    
    autodiff(ReverseWithPrimal, energy_sum!, Const(q), Duplicated(x, adpara.Fx), Duplicated(y, adpara.Fy), Duplicated(z, adpara.Fz), Const(interaction.n_atoms), Const(interaction.L), Const(interaction.α), Const(soepara), Duplicated(iterpara, adpara.iterpara_t), Const(interaction.k_set), Duplicated(adpara.U, adpara.dU))

    return nothing
end

function SoEwald2D_Fl!(interaction::SoEwald2DLongInteraction{T}, sys::MDSys, info::SimulationInfo{T}) where{T<:Number}

    revise_interaction!(interaction, sys, info)
    revise_adpara!(interaction.adpara, interaction.n_atoms)
    
    mass = interaction.mass
    ad_energy_sum!(interaction)

    for i in 1:sys.n_atoms
        info.particle_info[i].acceleration -= Point(interaction.adpara.Fx[i], interaction.adpara.Fy[i], interaction.adpara.Fz[i]) / (mass[i])
    end

    return nothing 
end