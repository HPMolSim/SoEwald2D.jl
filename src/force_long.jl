export ad_direct_sum, ad_energy_sum

function ad_direct_sum!(q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, TI}, adpara::AdPara{T}) where{T<:Number, TI<:Integer}
    revise_adpara!(adpara, para.n_atoms)
    autodiff(ReverseWithPrimal, direct_sum!, Const(q), Duplicated(x, adpara.Fx), Duplicated(y, adpara.Fy), Duplicated(z, adpara.Fz), Const(para), Duplicated(adpara.U, adpara.dU))
    return nothing
end