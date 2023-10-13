function update_iterpara_z!(iterpara::IterPara, z::Array{T}) where{T <: Number}
    sortperm!(iterpara.z_list, z)
    return nothing
end

# this function is used to update the z_list and m_list
function update_iterpara_m!(iterpara::IterPara, z::Array{T}, d::T) where{T <: Number}
    # d = k/2α²
    n_atoms = length(z)
    # here we manually define m_i = 1 if there are no particle out of the neighbor list
    # the out of neighbor members are 1:mᵢ - 1 and the neighbors are mᵢ:i - 1
    # for j in 1:m_i - 1
    # end
    # for j in m_i:i
    # end
    iterpara.m_list[1] = 1
    for i in 2:n_atoms
        m_i = iterpara.m_list[i - 1]
        while z[iterpara.z_list[i]] - z[iterpara.z_list[m_i]] > d && m_i < i
            m_i += 1
        end
        iterpara.m_list[i] = m_i
    end
    return nothing
end

function update_iterpara_A!(iterpara::IterPara, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, n_atoms::Int64, α::T, s::ComplexF64, K::Tuple{T, T, T}) where{T<:Number}
    k_x, k_y, k = K

    #update A
    iterpara.A[1] = zero(ComplexF64)
    for i in 2:n_atoms
        # qj exp (−ik · ρj − k zj + sl α zj)
        j = iterpara.z_list[i - 1]
        iterpara.A[i] = iterpara.A[i - 1] + q[j] * exp(- 1.0im * (k_x * x[j] + k_y * y[j]) + s * α * z[j])
    end

    return nothing
end

# notice that B have nothing to do with (s, w), so that it only needed to be updated once for each K
function update_iterpara_B!(iterpara::IterPara, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, n_atoms::Int64, α::T, K::Tuple{T, T, T}) where{T<:Number}
    k_x, k_y, k = K

    #update B
    iterpara.B[1] = zero(ComplexF64)
    for i in 2:n_atoms
        iterpara.B[i] = iterpara.B[i - 1]
        for m in iterpara.m_list[i - 1] : iterpara.m_list[i] - 1
            j = iterpara.z_list[m]
            iterpara.B[i] += 2 * q[j] * exp( - 1.0im * (k_x * x[j] + k_y * y[j]) + k * z[j])
        end 
    end

    return nothing
end

function update_iterpara_C!(iterpara::IterPara, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, n_atoms::Int64, α::T, s::ComplexF64, K::Tuple{T, T, T}) where{T<:Number}
    k_x, k_y, k = K

    #update C
    iterpara.C[1] = zero(ComplexF64)
    for i in 2:n_atoms
        iterpara.C[i] = iterpara.C[i - 1]
        for m in iterpara.m_list[i - 1] : iterpara.m_list[i] - 1
            j = iterpara.z_list[m]
            iterpara.C[i] += q[j] * exp( - 1.0im * (k_x * x[j] + k_y * y[j]) + k * z[j] + s * α * z[j])
        end 
    end

    return nothing
end

function update_iterpara_D!(iterpara::IterPara, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, n_atoms::Int64, α::T, s::ComplexF64, K::Tuple{T, T, T}) where{T<:Number}
    k_x, k_y, k = K
    n_atoms = n_atoms
    z_list = iterpara.z_list
    m_list = iterpara.m_list

    # n is the max value for the summation
    # n = para.α * para.L[3] * s
    if real(k - s * α) < 0
        # use the reverse mode iter
        iterpara.D[n_atoms] = zero(ComplexF64)
        for i in n_atoms - 1:-1:m_list[n_atoms]
            j = z_list[i]
            iterpara.D[n_atoms] += q[j] * exp( - 1.0im * (k_x * x[j] + k_y * y[j]) + k * z[j] - s * α * z[j])
        end

        for i in n_atoms - 1:-1:1
            li = z_list[i]
            iterpara.D[i] = iterpara.D[i + 1] - q[li] * exp(-1.0im * (k_x * x[li] + k_y * y[li]) + k * z[li] - s * α * z[li])
            for j in m_list[i + 1] - 1 : -1 : m_list[i]
                lj = z_list[j]
                iterpara.D[i] += q[lj] * exp(-1.0im * (k_x * x[lj] + k_y * y[lj]) + k * z[lj] - s * α * z[lj])
            end
        end
    else
        # use the forward mode iter
        iterpara.D[1] = zero(ComplexF64)
        for i in 2:n_atoms
            iterpara.D[i] = iterpara.D[i - 1]
            for j in m_list[i - 1]:m_list[i] - 1
                lj = z_list[j]
                iterpara.D[i] -= q[lj] * exp(-1.0im * (k_x * x[lj] + k_y * y[lj]) + k * z[lj] - s * α * z[lj])
            end
            li = z_list[i - 1]
            iterpara.D[i] += q[li] * exp(-1.0im * (k_x * x[li] + k_y * y[li]) + k * z[li] - s * α * z[li])
        end
    end
    
    return nothing
end

function energy_sum_k(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, n_atoms::Int64, α::T, soepara::SoePara{ComplexF64}, iterpara::IterPara) where{T<:Number}
    k_x, k_y, k = K
    update_iterpara_m!(iterpara, z, k / (2 * α^2))
    
    # compute the part irrelevant to SOE
    update_iterpara_B!(iterpara, q, x, y, z, n_atoms, α, K)
    sum_k = zero(ComplexF64)
    for i in 1:n_atoms
        sum_k += q[i] * q[i] * erfc(k/(2 * α))
        j = iterpara.z_list[i]
        sum_k += q[j] * exp(1.0im * (k_x * x[j] + k_y * y[j]) - k * z[j]) * iterpara.B[i]
    end

    #compute the SOE part
    for (s, w) in soepara.sw
        #update A, C and D
        update_iterpara_A!(iterpara, q, x, y, z, n_atoms, α, s, K)
        update_iterpara_C!(iterpara, q, x, y, z, n_atoms, α, s, K)
        update_iterpara_D!(iterpara, q, x, y, z, n_atoms, α, s, K)

        sum_k_soe_1 = zero(ComplexF64)
        sum_k_soe_2 = zero(ComplexF64)
        sum_k_soe_3 = zero(ComplexF64)
        
        for i in 1:n_atoms
            j = iterpara.z_list[i]
            sum_k_soe_1 += q[j] * exp(1.0im * (k_x * x[j] + k_y * y[j]) - s*α*z[j]) * iterpara.A[i]
            sum_k_soe_2 += - q[j] * exp(1.0im * (k_x * x[j] + k_y * y[j]) - k * z[j] - s*α*z[j]) * iterpara.C[i]
            sum_k_soe_3 += q[j] * exp(1.0im * (k_x * x[j] + k_y * y[j]) - k * z[j] + s*α*z[j]) * iterpara.D[i]
        end

        sum_k += s * w * exp(-k^2/(4 * α^2)) / (s * α + k) * α * sum_k_soe_1
        sum_k += w * (exp(s * k / (2 * α)) * sum_k_soe_2 + exp(- s * k / (2 * α)) * sum_k_soe_3)
    end
    return 2/k * real(sum_k)
end

function energy_sum_k0(q::Array{T}, z::Array{T}, n_atoms::Int64, α::T, soepara::SoePara{ComplexF64}, iterpara::IterPara) where{T<:Number}
    sum_k0 = zero(ComplexF64)
    z_list = iterpara.z_list

    for i in 1:n_atoms
        sum_k0 += q[i] * q[i] * (1 / (α * sqrt(π)))
    end

    iterpara.A[1] = zero(ComplexF64)
    iterpara.B[1] = zero(ComplexF64)
    iterpara.C[1] = zero(ComplexF64)
    iterpara.D[1] = zero(ComplexF64)

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

    # relevant to (s, w)
    for (s, w) in soepara.sw
        soe_sum_k0_1 = zero(ComplexF64) # for the erfc part
        soe_sum_k0_2 = zero(ComplexF64) # for the exp(-(α z_ij)^2) part
        for i in 2:n_atoms
            l = z_list[i - 1]
            t = q[l] * exp(s * α * z[l])
            iterpara.C[i] = iterpara.C[i - 1] + t
            iterpara.D[i] = iterpara.D[i - 1] + z[l] * t
        end
        for i in 1:n_atoms
            l = z_list[i]
            q_i = q[l]
            z_i = z[l]
            soe_sum_k0_1 += q_i * (z_i * iterpara.C[i] - iterpara.D[i]) * exp( - s * α * z_i)
            soe_sum_k0_2 += q_i * iterpara.C[i] * exp( - s * α * z_i)
        end
        sum_k0 -= 2 * w * soe_sum_k0_1 
        sum_k0 += 2 * w * s / (α * sqrt(π)) * soe_sum_k0_2 * sqrt(π) / 2
    end
    return real(sum_k0)
end

function energy_sum!(q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, n_atoms::Int64, L::NTuple{3, T}, α::T, soepara::SoePara, iterpara::IterPara, k_set::Array{NTuple{3, T}}, U::Array{T}) where{T<:Number}
    energy = zero(T)

    update_iterpara_z!(iterpara, z)

    U_k0 = - energy_sum_k0(q, z, n_atoms, α, soepara, iterpara) * π / (L[1] * L[2])

    for i in 1:size(k_set, 1)
        K = k_set[i]
        energy += energy_sum_k(K, q, x, y, z, n_atoms, α, soepara, iterpara)
    end

    energy *= π / (2 * L[1] * L[2])
    energy += U_k0
    U[1] = energy / 4π
    return nothing
end

# # the summation will be summed up first, and the a interface SoEwald2D_EL will be added
function SoEwald2D_El(interaction::SoEwald2DLongInteraction{T}, sys::MDSys, info::SimulationInfo{T}) where{T<:Number}
    
    U = [zero(T)]

    revise_interaction!(interaction, sys, info)
    energy_sum!(interaction.q, interaction.x, interaction.y, interaction.z, interaction.n_atoms, interaction.L, interaction.α, interaction.soepara, interaction.iterpara, interaction.k_set, U)
    
    return U[1]
end