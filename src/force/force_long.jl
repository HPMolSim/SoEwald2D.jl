Base.real(x::Point) = Point(real.(x.coo))

function energy_sum_k!(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, n_atoms::Int64, α::T, soepara::SoePara{ComplexF64}, iterpara::IterPara, U::Array{T}) where{T<:Number}
    U[1] = energy_sum_k(K, q, x, y, z, n_atoms, α, soepara, iterpara)
    return nothing
end

function energy_sum_k0!(q::Array{T}, z::Array{T}, n_atoms::Int64, α::T, soepara::SoePara{ComplexF64}, iterpara::IterPara, U::Array{T}) where{T<:Number}
    U[1] = energy_sum_k0(q, z, n_atoms, α, soepara, iterpara)
    return nothing
end

function force_sum_k(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, n_atoms::Int64, α::T, soepara::SoePara{ComplexF64}, iterpara::IterPara, adpara::AdPara) where{T<:Number}

    revise_adpara!(adpara, n_atoms)

    autodiff(ReverseWithPrimal, energy_sum_k!, Const(K), Const(q), Duplicated(x, adpara.Fx), Duplicated(y, adpara.Fy), Duplicated(z, adpara.Fz), Const(n_atoms), Const(α), Const(soepara), Duplicated(iterpara, adpara.iterpara_t), Duplicated(adpara.U, adpara.dU))

    return [adpara.Fx, adpara.Fy, adpara.Fz]
end

function force_sum_k0(q::Array{T}, z::Array{T}, n_atoms::Int64, α::T, soepara::SoePara{ComplexF64}, iterpara::IterPara, adpara::AdPara) where{T<:Number}

    revise_adpara!(adpara, n_atoms)
    autodiff(ReverseWithPrimal, energy_sum_k0!, Const(q), Duplicated(z, adpara.Fz), Const(n_atoms), Const(α), Const(soepara), Duplicated(iterpara, adpara.iterpara_t), Duplicated(adpara.U, adpara.dU))

    return [adpara.Fx, adpara.Fy, adpara.Fz]
end

function force_sum(interaction::SoEwald2DLongInteraction{T}) where{T}

    iterpara = interaction.iterpara
    soepara = interaction.soepara
    adpara = interaction.adpara
    q = interaction.q
    x = interaction.x
    y = interaction.y
    z = interaction.z

    rbm = interaction.rbm
    rbm_p = interaction.rbm_p
    P = interaction.P

    α = interaction.α
    n_atoms = interaction.n_atoms
    L = interaction.L
    k_set = interaction.k_set
    ϵ_0 = interaction.ϵ_0
    parallel = interaction.parallel
    rng = interaction.rng
    
    update_iterpara_z!(iterpara, z)

    F_k0 = - force_sum_k0(q, z, n_atoms, α, soepara, iterpara, adpara) * π / (L[1] * L[2])
    F_k = [zeros(T, n_atoms), zeros(T, n_atoms), zeros(T, n_atoms)]

    if parallel
        if rbm == false
            F_k = @distributed (+) for k in k_set
                exp(- k[3]^2 / (4 * α^2)) * force_sum_k(k, q, x, y, z, n_atoms, α, soepara, iterpara, adpara)
            end
        else
            F_k = @distributed (+) for i in 1:rbm_p
                P / rbm_p * force_sum_k(k_set[rand(rng, 1:end)], q, x, y, z, n_atoms, α, soepara, iterpara, adpara)
            end
        end
    else
        if rbm == false
            for k in k_set
                F_k += exp(- k[3]^2 / (4 * α^2)) * force_sum_k(k, q, x, y, z, n_atoms, α, soepara, iterpara, adpara)
            end
        else
            for i in 1:rbm_p
                F_k += P / rbm_p * force_sum_k(k_set[rng, rand(1:end)], q, x, y, z, n_atoms, α, soepara, iterpara, adpara)
            end
        end
    end

    F_k *= π / (2 * L[1] * L[2])
    F_k .+= F_k0

    return F_k / (4π * ϵ_0)
end

function SoEwald2D_Fl!(interaction::SoEwald2DLongInteraction{T}, sys::MDSys, info::SimulationInfo{T}) where{T<:Number}

    revise_interaction!(interaction, sys, info)
    revise_adpara!(interaction.adpara, interaction.n_atoms)
    
    mass = interaction.mass
    Fx, Fy, Fz = force_sum(interaction)

    for i in 1:sys.n_atoms
        info.particle_info[i].acceleration -= Point(Fx[i], Fy[i], Fz[i]) / (mass[i])
    end

    return nothing 
end