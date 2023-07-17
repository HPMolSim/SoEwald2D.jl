function direct_sum_S1(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, TI}) where{T<:Number, TI<:Integer}
    k_x, k_y, k = K
    sum_k = big(zero(T))
    α = para.α
    for i in 1:para.n_atoms
        for j in 1:i - 1
            x_ij = x[i] - x[j]
            y_ij = y[i] - y[j]
            z_ij = abs(z[i] - z[j])
            sum_k += q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * exp(big(k * z_ij)) * erfc(big(k / (2α) + α * z_ij))
        end
    end
    return 2 * sum_k / k
end

function direct_sum_S2(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, TI}) where{T<:Number, TI<:Integer}
    k_x, k_y, k = K
    sum_k0 = zero(T)
    sum_k1 = zero(T)
    sum_k2 = zero(T)
    α = para.α
    for i in 1:para.n_atoms
        for j in 1:para.n_atoms
            if i ≠ j
                x_ij = x[i] - x[j]
                y_ij = y[i] - y[j]
                z_ij =  - abs(z[i] - z[j])
                if k / (2α) + α * z_ij > 0
                    sum_k2 += q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * (exp(k * z_ij) * erfc(k / (2α) + α * z_ij)) # correspond to D
                else
                    sum_k0 += q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * exp(k * z_ij) * 2 # correspond to B
                    sum_k1 += - q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * (exp(k * z_ij) * erfc( - k / (2α) -  α * z_ij)) # correspond to C
                end
            end
        end
    end
    return sum_k0 / k, sum_k1 / k, sum_k2 / k
end

function direct_sum_B(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, TI}) where{T<:Number, TI<:Integer}
    n_atoms = para.n_atoms
    k_x, k_y, k = K
    z_list = sortperm(z)
    d = k / (2 * para.α^2)
    B = zeros(ComplexF64, n_atoms)

    for i in 2:n_atoms
        li = z_list[i]
        for j in 1:i - 1
            lj = z_list[j]
            dz = z[li] - z[lj]
            if dz > d
                B[i] += 2 * q[lj] * exp( - 1.0im * (k_x * x[lj] + k_y * y[lj]) + k * z[lj])
            else
                break
            end
        end 
    end
    return B
end

function direct_sum_CD(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, TI}, s::ComplexF64) where{T<:Number, TI<:Integer}
    n_atoms = para.n_atoms
    k_x, k_y, k = K
    z_list = sortperm(z)
    d = k / (2 * para.α^2)
    C = zeros(ComplexF64, n_atoms)
    D = zeros(ComplexF64, n_atoms)

    for i in 2:n_atoms
        li = z_list[i]
        for j in 1:i - 1
            lj = z_list[j]
            dz = z[li] - z[lj]
            if dz > d
                C[i] += q[lj] * exp( - 1.0im * (k_x * x[lj] + k_y * y[lj]) + (k + s * para.α) * z[lj])
            else
                D[i] += q[lj] * exp( - 1.0im * (k_x * x[lj] + k_y * y[lj]) + (k - s * para.α) * z[lj])
            end
        end 
    end
    return C, D
end

