# Validation reference for the long-range FORCE. Despite the name, what this
# returns is `-∇U` -- the same sign as `SoEwald2D.force` -- not the gradient:
# the pre-decoupling test compared it directly against
# `info.particle_info[i].acceleration` after `SoEwald2D_Fl!`, which is
# `-∇U / mass` at unit mass. Same units convention as direct_sum.jl: no
# 1 / (4π ϵ_0) prefactor.
function diff_direct_sum(plan::SoEwald2DLongPlan{T}, poses, charges) where{T<:Number}

    n_atoms = plan.n_atoms
    
    sum = [SVector{3, T}(zero(T), zero(T), zero(T)) for _=1:n_atoms]

    _scatter_long!(plan, poses, charges)
    diff_direct_sum_k0!(plan.q, plan.z, plan, sum)
    for K in plan.k_set
        diff_direct_sum_k!(K, plan.q, plan.x, plan.y, plan.z, plan, sum)
    end
    return - sum .* T(2)
end

function diff_direct_sum_k0!(q::Array{T}, z::Array{T}, para::SoEwald2DLongPlan{T}, sum::Vector{SVector{3, T}}) where {T<:Number}
    α = para.α
    for i in 1:para.n_atoms
        for j in 1:para.n_atoms
            z_ij = z[i] - z[j]
            sum[i] -= SVector{3, T}(zero(T), zero(T), q[i] * q[j] * (erf(α * z_ij))  / (4 *  para.L[1] * para.L[2]))
        end
    end
    return nothing
end

function diff_direct_sum_k!(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DLongPlan{T},sum::Vector{SVector{3, T}}) where{T<:Number}

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
        sum[i] += SVector{3, T}(sum_x, sum_y, sum_z) / (8 * para.L[1] * para.L[2])
    end

    return nothing
end