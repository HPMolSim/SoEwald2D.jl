function diff_direct_sum(interaction::SoEwald2DLongInteraction{T}, sys::MDSys{T}, info::SimulationInfo{T}) where{T<:Number}

    n_atoms = interaction.n_atoms
    
    sum = [Point(zero(T), zero(T), zero(T)) for _=1:n_atoms]

    revise_interaction!(interaction, sys, info)
    diff_direct_sum_k0!(interaction.q, interaction.z, interaction, sum)
    for K in interaction.k_set
        diff_direct_sum_k!(K, interaction.q, interaction.x, interaction.y, interaction.z, interaction, sum)
    end
    return - sum .* T(2)
end

function diff_direct_sum_k0!(q::Array{T}, z::Array{T}, para::SoEwald2DLongInteraction{T}, sum::Vector{Point{3, T}}) where {T<:Number}
    α = para.α
    for i in 1:para.n_atoms
        for j in 1:para.n_atoms
            z_ij = z[i] - z[j]
            sum[i] -= Point(zero(T), zero(T), q[i] * q[j] * (erf(α * z_ij))  / (4 *  para.L[1] * para.L[2]))
        end
    end
    return nothing
end

function diff_direct_sum_k!(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DLongInteraction{T},sum::Vector{Point{3, T}}) where{T<:Number}

    n_atoms = para.n_atoms

    k_x, k_y, k = K
    α = para.α
    for i in 1:para.n_atoms
        sum_x = zero(T)
        sum_y = zero(T)
        sum_z = zero(T)
        for j in 1:para.n_atoms
            x_ij = x[i] - x[j]
            y_ij = y[i] - y[j]
            z_ij = z[i] - z[j]

            sum_xy = - q[i] * q[j] * sin(k_x * x_ij + k_y * y_ij) * (exp(k * z_ij) * erfc(k / (2α) + α * z_ij) + exp( - k * z_ij) * erfc(k / (2α) - α * z_ij))
            sum_x += k_x * sum_xy / k
            sum_y += k_y * sum_xy / k

            sum_z += q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * (
                k * exp(k * z_ij) * erfc(k / (2α) + α * z_ij) - 
                k * exp( - k * z_ij) * erfc(k / (2α) - α * z_ij) -
                2α / sqrt(π) * exp(k * z_ij) * exp(-(k / (2α) + α * z_ij)^2) +
                2α / sqrt(π) * exp(- k * z_ij) * exp(-(k / (2α) - α * z_ij)^2) ) / k
        end
        sum[i] += Point(sum_x, sum_y, sum_z) / (8 * para.L[1] * para.L[2])
    end

    return nothing
end