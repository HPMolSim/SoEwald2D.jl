export ad_energy_sum!

function ad_energy_sum!(q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, TI}, soepara::SoePara{ComplexF64}, iterpara::IterPara, adpara::AdPara{T}) where{T<:Number, TI<:Integer}
    revise_adpara!(adpara, para.n_atoms)
    autodiff(ReverseWithPrimal, energy_sum!, Const(q), Duplicated(x, adpara.Fx), Duplicated(y, adpara.Fy), Duplicated(z, adpara.Fz), Const(para), Const(soepara), Duplicated(iterpara, adpara.iterpara_t), Duplicated(adpara.U, adpara.dU))

    return nothing
end