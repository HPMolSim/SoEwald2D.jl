export direct_sum!, direct_sum_k, diff_direct_sum!, ad_direct_sum!, direct_sum_k_ABCD_soerfc, direct_sum_k_ABCD, direct_iterpara_D, direct_iter_sum_D

# this function will calculate the direct sum
function direct_sum!(q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, TI}, U::Vector{T}) where{T<:Number, TI<:Integer}
    U[1] = zero(T)

    U_k0 = - direct_sum_k0(q, z, para) * π / (para.L[1] * para.L[2])

    for K in para.k_set
        U[1] += direct_sum_k(K, q, x, y, z, para)
    end

    U[1] *= π / (2 * para.L[1] * para.L[2])
    U[1] += U_k0
    return nothing
end

function direct_sum_k(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, TI}) where{T<:Number, TI<:Integer}
    k_x, k_y, k = K
    sum_k = zero(T)
    α = para.α
    for i in 1:para.n_atoms
        for j in 1:para.n_atoms
            x_ij = x[i] - x[j]
            y_ij = y[i] - y[j]
            z_ij = z[i] - z[j]
            sum_k += q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * (exp(k * z_ij) * erfc(k / (2α) + α * z_ij) + exp( - k * z_ij) * erfc(k / (2α) - α * z_ij))
        end
    end
    return sum_k / k
end

function direct_sum_k_ABCD(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, TI}) where{T<:Number, TI<:Integer}
    k_x, k_y, k = K
    sum_k_A = zero(T)
    sum_k_B = zero(T)
    sum_k_C = zero(T)
    sum_k_D = zero(T)
    α = para.α
    for i in 1:para.n_atoms
        for j in 1:para.n_atoms
            if i ≠ j
                x_ij = x[i] - x[j]
                y_ij = y[i] - y[j]
                z_ij = z[i] - z[j]

                sum_k_A += q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * exp(k * abs(z_ij)) * erfc(k / (2α) + α * abs(z_ij))
                if k / (2α) - α * abs(z_ij) < 0
                    sum_k_B += 2 * q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * exp( - k * abs(z_ij))
                    sum_k_C += - q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * exp( - k * abs(z_ij)) * erfc( - k / (2α) + α * abs(z_ij))
                else
                    sum_k_D += q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * exp( - k * abs(z_ij)) * erfc(k / (2α) - α * abs(z_ij))
                end
            end
        end
    end
    return sum_k_A / k, sum_k_B / k, sum_k_C / k, sum_k_D / k
end

function direct_sum_k_ABCD_soerfc(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, TI}, soepara::SoePara{ComplexF64}) where{T<:Number, TI<:Integer}
    k_x, k_y, k = K
    sum_k_A = zero(T)
    sum_k_B = zero(T)
    sum_k_C = zero(T)
    sum_k_D = zero(T)
    α = para.α
    for i in 1:para.n_atoms
        for j in 1:para.n_atoms
            if i ≠ j
                x_ij = x[i] - x[j]
                y_ij = y[i] - y[j]
                z_ij = z[i] - z[j]

                sum_k_A += q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * exp(k * abs(z_ij)) * soerfc(k / (2α) + α * abs(z_ij), soepara)
                if k / (2α) - α * abs(z_ij) < 0
                    sum_k_B += 2 * q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * exp( - k * abs(z_ij))
                    sum_k_C += - q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * exp( - k * abs(z_ij)) * soerfc( - k / (2α) + α * abs(z_ij), soepara)
                else
                    sum_k_D += q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * exp( - k * abs(z_ij)) * soerfc(k / (2α) - α * abs(z_ij), soepara)
                end
            end
        end
    end
    return sum_k_A / k, sum_k_B / k, sum_k_C / k, sum_k_D / k
end

function direct_iterpara_D(q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, Int64}, s::ComplexF64, K::Tuple{T, T, T}) where{T<:Number}
    k_x, k_y, k = K
    n_atoms = para.n_atoms
    α = para.α

    D = zeros(ComplexF64, n_atoms)
    z_list = sortperm(z)

    for i in 1:n_atoms
        li = z_list[i]
        for j in 1: i - 1
            lj = z_list[j]
            z_ij = z[li] - z[lj]
            if k / (2α) - α * abs(z_ij) > 0
                D[i] += q[lj] * exp( - 1.0im * (k_x * x[lj] + k_y * y[lj]) + k * z[lj] - s * para.α * z[lj])
            end
        end
    end

    return D
end

function direct_iter_sum_D(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, Int64}, soepara::SoePara{ComplexF64}, iterpara::IterPara) where{T<:Number}
    k_x, k_y, k = K
    n_atoms = para.n_atoms
    α = para.α
    update_iterpara_z!(iterpara, z)
    update_iterpara_m!(iterpara, z, k / (2 * α^2))

    sum_k1 = zero(ComplexF64)
    sum_k2 = zero(ComplexF64)

    for (s, w) in soepara.sw
        D_dir = direct_iterpara_D(q, x, y, z, para, s, K)
        update_iterpara_D!(iterpara, q, x, y, z, para, s, K)

        sum_k_soe_1 = zero(ComplexF64)
        sum_k_soe_2 = zero(ComplexF64)
        
        for i in 1:para.n_atoms
            j = iterpara.z_list[i]
            sum_k_soe_1 += q[j] * exp( + 1.0im * (k_x * x[j] + k_y * y[j])) * exp(- k * z[j] + s*para.α*z[j] - s * k / (2 * para.α)) * iterpara.D[i]
            sum_k_soe_2 += q[j] * exp( + 1.0im * (k_x * x[j] + k_y * y[j])) * exp(- k * z[j] + s*para.α*z[j] - s * k / (2 * para.α)) * D_dir[i]
            if abs(sum_k_soe_1 - sum_k_soe_2) > 1e-10
                temp = q[j] * exp( + 1.0im * (k_x * x[j] + k_y * y[j])) * exp(- k * z[j] + s*para.α*z[j] - s * k / (2 * para.α))
                @show iterpara.D[i], D_dir[i], sum_k_soe_1, sum_k_soe_2, temp
            end
        end

        sum_k1 += w * sum_k_soe_1
        sum_k2 += w * sum_k_soe_2
    end
    return 2/k * real(sum_k1), 2/k * real(sum_k2)
end

function direct_sum_k0(q::Array{T}, z::Array{T}, para::SoEwald2DPara{T, TI}) where{T<:Number, TI<:Integer}
    sum_k0 = zero(T)
    α = para.α
    for i in 1:para.n_atoms
        for j in 1:para.n_atoms
            z_ij = z[i] - z[j]
            sum_k0 += q[i] * q[j] * (1 / (α * sqrt(π)) * exp(-(α * z_ij)^2) + z_ij * erf(α * z_ij))
        end
    end
    return sum_k0
end

function diff_direct_sum!(q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, TI}, F_x::Vector{T}, F_y::Vector{T}, F_z::Vector{T}) where{T<:Number, TI<:Integer}

    diff_direct_sum_k0!(q, z, para, F_z)
    for K in para.k_set
        diff_direct_sum_k!(K, q, x, y, z, para, F_x, F_y, F_z)
    end
    return nothing
end

function diff_direct_sum_k0!(q::Array{T}, z::Array{T}, para::SoEwald2DPara{T, TI}, F_z::Vector{T}) where {T<:Number, TI<:Integer}
    α = para.α
    for i in 1:para.n_atoms
        for j in 1:para.n_atoms
            z_ij = z[i] - z[j]
            F_z[i] -= q[i] * q[j] * (erf(α * z_ij)) * π / (para.L[1] * para.L[2])
        end
    end
    return nothing
end

function diff_direct_sum_k!(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, TI}, F_x::Vector{T}, F_y::Vector{T}, F_z::Vector{T}) where{T<:Number, TI<:Integer}
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
        F_x[i] += sum_x * π / (2 * para.L[1] * para.L[2])
        F_y[i] += sum_y * π / (2 * para.L[1] * para.L[2])
        F_z[i] += sum_z * π / (2 * para.L[1] * para.L[2])
    end

    return nothing
end

function ad_direct_sum!(q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, TI}, adpara::AdPara{T}) where{T<:Number, TI<:Integer}
    revise_adpara!(adpara, para.n_atoms)
    autodiff(ReverseWithPrimal, direct_sum!, Const(q), Duplicated(x, adpara.Fx), Duplicated(y, adpara.Fy), Duplicated(z, adpara.Fz), Const(para), Duplicated(adpara.U, adpara.dU))
    return nothing
end


function f1(z, α)
    return 1/(α * sqrt(π)) * exp(-(α * z)^2) + z * erf(α * z)
end

function diff_f1(z, α)
    return erf(α * z)
end