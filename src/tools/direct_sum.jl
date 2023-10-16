# this function will calculate the direct sum
function direct_sum(interaction::SoEwald2DLongInteraction{T}, sys::MDSys{T}, info::SimulationInfo{T}) where{T<:Number}

    energy = zero(T)
    for i in 1:interaction.n_atoms
        energy += direct_sum(i, interaction, sys, info)
    end
    return energy
end

function direct_sum(i::Int, interaction::SoEwald2DLongInteraction{T}, sys::MDSys{T}, info::SimulationInfo{T}) where{T<:Number}

    energy = zero(T)
    revise_interaction!(interaction, sys, info)

    q = interaction.q
    x = interaction.x
    y = interaction.y
    z = interaction.z

    U_k0 = - direct_sum_k0i(i, q, z, interaction)

    for K in interaction.k_set
        energy += direct_sum_ki(i, K, q, x, y, z, interaction)
    end

    energy += U_k0
    return energy
end


function direct_sum_ki(i::Int, K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DLongInteraction{T}) where{T<:Number}
    k_x, k_y, k = K
    α = para.α
    sum_ki = zero(T)
    for j in 1:para.n_atoms
        x_ij = x[i] - x[j]
        y_ij = y[i] - y[j]
        z_ij = z[i] - z[j]
        sum_ki += q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * (exp(k * z_ij) * erfc(k / (2α) + α * z_ij) + exp( - k * z_ij) * erfc(k / (2α) - α * z_ij)) / (8 * para.L[1] * para.L[2] * k)
    end
    return sum_ki
end


function direct_sum_k0i(i::Int, q::Array{T}, z::Array{T}, para::SoEwald2DLongInteraction{T}) where{T<:Number}
    α = para.α
    sum_k0i = zero(T)
    for j in 1:para.n_atoms
        z_ij = z[i] - z[j]
        sum_k0i += q[i] * q[j] * (1 / (α * sqrt(π)) * exp(-(α * z_ij)^2) + z_ij * erf(α * z_ij)) / (4 *  para.L[1] * para.L[2])
    end
    return sum_k0i
end

function soe_direct_sum(interaction::SoEwald2DLongInteraction{T}, sys::MDSys{T}, info::SimulationInfo{T}, soepara::SoePara{ComplexF64}) where{T<:Number}

    energy = zero(T)
    for i in 1:interaction.n_atoms
        energy += soe_direct_sum(i, interaction, sys, info, soepara)
    end
    return energy
end

function soe_direct_sum(i::Int, interaction::SoEwald2DLongInteraction{T}, sys::MDSys{T}, info::SimulationInfo{T}, soepara::SoePara{ComplexF64}) where{T<:Number}

    energy = zero(T)
    revise_interaction!(interaction, sys, info)

    q = interaction.q
    x = interaction.x
    y = interaction.y
    z = interaction.z

    U_k0 = - soe_direct_sum_k0i(i, q, z, interaction, soepara)

    for K in interaction.k_set
        energy += soe_direct_sum_ki(i, K, q, x, y, z, interaction, soepara)
    end

    energy += U_k0
    return energy
end

function soe_direct_sum_ki(i::Int, K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DLongInteraction{T}, soepara::SoePara{ComplexF64}) where{T<:Number}
    k_x, k_y, k = K
    α = para.α
    sum_ki = zero(ComplexF64)
    for j in 1:para.n_atoms
        x_ij = x[i] - x[j]
        y_ij = y[i] - y[j]
        z_ij = z[i] - z[j]
        z_ij = abs(z_ij)
        sum_ki += q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * (soexp_mul_erfc(z_ij, α, k, soepara) + soexp_mul_erfc(-z_ij, α, k, soepara)) / (8 * para.L[1] * para.L[2] * k)
    end
    return sum_ki
end

function soe_direct_sum_k0i(i::Int, q::Array{T}, z::Array{T}, para::SoEwald2DLongInteraction{T}, soepara::SoePara{ComplexF64}) where{T<:Number}
    α = para.α
    sum_k0i = zero(ComplexF64)
    for j in 1:para.n_atoms
        z_ij = z[i] - z[j]
        sum_k0i += q[i] * q[j] * (1 / (α * sqrt(π)) * soexp(-α * z_ij, soepara) + z_ij * soerf(α * z_ij, soepara)) / (4 *  para.L[1] * para.L[2])
    end
    return sum_k0i
end