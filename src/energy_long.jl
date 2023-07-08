export SoEwald2D_El, energy_sum, direct_sum!

# this function will calculate the direct sum
# the summation is given by
# 
function direct_sum!(q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, TI}, U::Vector{T}) where{T<:Number, TI<:Integer}
    U[1] = zero(T)
    α = para.α
    for (k_x, k_y, k) in para.k_set
        sum_k = zero(T)
        for i in 1:para.n_atoms
            for j in 1:para.n_atoms
                x_ij = x[i] - x[j]
                y_ij = y[i] - y[j]
                z_ij = z[i] - z[j]
                sum_k += q[i] * q[j] * cos(k_x * x_ij + k_y * y_ij) * (exp(k * z_ij) * erfc(k / (2α) + α * z_ij) + exp( - k * z_ij) * erfc(k / (2α) - α * z_ij))
            end
        end
        U[1] += sum_k / k
    end
    U[1] *= π / (2 * para.L[1] * para.L[2])
    return nothing
end

# # the summation will be summed up first, and the a interface SoEwald2D_EL will be added
function energy_sum!(q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2D{T, TI}, soepara::SoePara{TC}, U::Vector{T}) where{T<:Number, TI<:Integer, TC<:Complex}
    U[1] = zero(T)
    
end