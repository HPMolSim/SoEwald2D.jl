function update_iterpara_z!(iterpara::IterPara, z::Array{T}) where{T <: Number}
    sortperm!(iterpara.z_list, z)
    return nothing
end

function update_iterpara_A!(iterpara::IterPara, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, n_atoms::Int64, α::T, s::ComplexF64, K::Tuple{T, T, T}) where{T<:Number}
    k_x, k_y, k = K

    #update A
    iterpara.A[1] = zero(ComplexF64)
    j0 = iterpara.z_list[1]
    iterpara.A[2] = q[j0] * exp(- 1.0im * (k_x * x[j0] + k_y * y[j0]))
    for i in 3:n_atoms
        # qj exp (−ik · ρj − k zj + sl α zj)
        j = iterpara.z_list[i - 1]
        l = iterpara.z_list[i - 2]
        iterpara.A[i] = iterpara.A[i - 1] * exp(s * α * (z[l] - z[j])) + q[j] * exp(- 1.0im * (k_x * x[j] + k_y * y[j]))
    end

    return nothing
end

# notice that B have nothing to do with (s, w), so that it only needed to be updated once for each K
function update_iterpara_B!(iterpara::IterPara, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, n_atoms::Int64, α::T, K::Tuple{T, T, T}) where{T<:Number}
    k_x, k_y, k = K

    #update B
    j0 = iterpara.z_list[1]
    iterpara.B[2] = q[j0] * exp(- 1.0im * (k_x * x[j0] + k_y * y[j0]))

    for i in 3:n_atoms
        j = iterpara.z_list[i - 1]
        l = iterpara.z_list[i - 2]
        iterpara.B[i] = iterpara.B[i - 1] * exp(k * (z[l] - z[j])) + q[j] * exp(- 1.0im * (k_x * x[j] + k_y * y[j]))
    end

    return nothing
end

function energy_sum_k(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, n_atoms::Int64, α::T, soepara::SoePara{ComplexF64}, iterpara::IterPara) where{T<:Number}
    k_x, k_y, k = K

    # constants for the part irrelevant to SOE
    c_1 = zero(ComplexF64)
    c_2 = zero(ComplexF64)
    for (s, w) in soepara.sw
        c_1 += α * w * 2 / sqrt(π) / (s * α + k)
        c_2 += 2.0 * α * w * s * 2.0 / sqrt(π) * α / ((s * α)^2 - k^2)
    end
    
    # compute the part irrelevant to SOE
    update_iterpara_B!(iterpara, q, x, y, z, n_atoms, α, K)
    sum_k = zero(ComplexF64)
    for i in 1:n_atoms
        sum_k += q[i] * q[i] * c_1
    end

    for i in 2:n_atoms
        j = iterpara.z_list[i]
        l = iterpara.z_list[i - 1]
        sum_k += c_2 * q[j] * exp(1.0im * (k_x * x[j] + k_y * y[j]) - k * (z[j] - z[l])) * iterpara.B[i]
    end

    #compute the SOE part
    for (s, w) in soepara.sw
        update_iterpara_A!(iterpara, q, x, y, z, n_atoms, α, s, K)

        sum_k_soe = zero(ComplexF64)
        
        for i in 2:n_atoms
            j = iterpara.z_list[i]
            l = iterpara.z_list[i - 1]
            sum_k_soe += q[j] * exp(1.0im * (k_x * x[j] + k_y * y[j]) - s*α*z[j] + s*α*z[l]) * iterpara.A[i]
        end

        sum_k += 2.0 * k * 2.0 / sqrt(π) * w * α * sum_k_soe / (k^2.0 - (s * α)^2.0)
    end

    return 2.0 / k * real(sum_k)
end

function energy_sum_k0(q::Array{T}, z::Array{T}, n_atoms::Int64, α::T, soepara::SoePara{ComplexF64}, iterpara::IterPara) where{T<:Number}
    sum_k0 = zero(ComplexF64)
    z_list = iterpara.z_list

    for i in 1:n_atoms
        sum_k0 += q[i] * q[i] * (1 / (α * sqrt(π)))
    end

    iterpara.A[1] = zero(ComplexF64)
    iterpara.B[1] = zero(ComplexF64)

    # irrelevant to (s, w)
    for i in 2:n_atoms
        l = z_list[i - 1]
        iterpara.A[i] = iterpara.A[i - 1] + q[l]
        iterpara.B[i] = iterpara.B[i - 1] + q[l] * z[l]
    end
    for i in 1:n_atoms
        l = z_list[i]
        q_i = q[l]
        z_i = z[l]
        sum_k0 += 2 * q_i * (z_i * iterpara.A[i] - iterpara.B[i])
    end

    iterpara.A[1] = zero(ComplexF64)
    iterpara.B[1] = zero(ComplexF64)
    # relevant to (s, w)
    for (s, w) in soepara.sw
        soe_sum_k0_1 = zero(ComplexF64) # for the erfc part
        soe_sum_k0_2 = zero(ComplexF64) # for the exp(-(α z_ij)^2) part
        j0 = z_list[1]
        iterpara.A[2] = q[j0]
        iterpara.B[2] = z[j0] * q[j0]
        for i in 3:n_atoms
            j = z_list[i - 1]
            l = z_list[i - 2]
            iterpara.A[i] = iterpara.A[i - 1] * exp(s * α * (z[l] - z[j])) + q[j]
            iterpara.B[i] = iterpara.B[i - 1] * exp(s * α * (z[l] - z[j])) + z[j] * q[j]
        end
        for i in 2:n_atoms
            j = z_list[i]
            l = z_list[i - 1]
            soe_sum_k0_1 += q[j] * (z[j] * iterpara.A[i] - iterpara.B[i]) * exp( - s * α * (z[j] - z[l]))
            soe_sum_k0_2 += q[j] * iterpara.A[i] * exp( - s * α * (z[j] - z[l]))
        end
        sum_k0 -= 2 * w / s * 2.0 / sqrt(π) * soe_sum_k0_1 
        sum_k0 += 2 * w / (α * sqrt(π)) * soe_sum_k0_2
    end
    return real(sum_k0)
end

function energy_sum!(q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, n_atoms::Int64, ϵ_0::T, L::NTuple{3, T}, α::T, soepara::SoePara, iterpara::IterPara, k_set::Array{NTuple{3, T}}, rbm::Bool, rbm_p::Int, P::T, U::Array{T}) where{T<:Number}
    energy = zero(T)

    update_iterpara_z!(iterpara, z)

    U_k0 = - energy_sum_k0(q, z, n_atoms, α, soepara, iterpara) * π / (L[1] * L[2])

    if rbm == false
        for i in 1:size(k_set, 1)
            K = k_set[i]
            k = K[3]
            energy += exp(- k^2 / (4 * α^2)) * energy_sum_k(K, q, x, y, z, n_atoms, α, soepara, iterpara)
        end
    else
        for i in 1:rbm_p
            K = k_set[rand(1:size(k_set, 1))]
            energy += P / rbm_p * energy_sum_k(K, q, x, y, z, n_atoms, α, soepara, iterpara)
        end
    end

    energy *= π / (2 * L[1] * L[2])
    energy += U_k0
    U[1] = energy / (4π * ϵ_0)
    return nothing
end

# # the summation will be summed up first, and the a interface SoEwald2D_EL will be added
function SoEwald2D_El(interaction::SoEwald2DLongInteraction{T}, sys::MDSys, info::SimulationInfo{T}) where{T<:Number}
    
    U = [zero(T)]

    revise_interaction!(interaction, sys, info)
    energy_sum!(interaction.q, interaction.x, interaction.y, interaction.z, interaction.n_atoms, interaction.ϵ_0, interaction.L, interaction.α, interaction.soepara, interaction.iterpara, interaction.k_set, interaction.rbm, interaction.rbm_p, interaction.P, U)
    
    return U[1]
end