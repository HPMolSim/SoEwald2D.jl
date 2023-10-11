Base.real(x::Point) = Point(real.(x.coo))

function update_iterpara_A12!(iterpara::IterPara, q::Vector{T}, x::Vector{T}, y::Vector{T}, z::Vector{T}, para::SoEwald2DLongInteraction{T}, s::ComplexF64, K::Tuple{T, T, T}) where{T<:Number}
    k_x, k_y, k = K

    #update A
    iterpara.A1[1] = zero(ComplexF64)
    iterpara.A2[para.n_atoms] = zero(ComplexF64)
    for i in 2:para.n_atoms
        # qj exp (−ik · ρj − k zj + sl α zj)
        j = iterpara.z_list[i - 1]
        iterpara.A1[i] = iterpara.A1[i - 1] + q[j] * exp(- 1.0im * (k_x * x[j] + k_y * y[j]) + s * para.α * z[j])
    end

    for i in para.n_atoms - 1:-1:1
        j = iterpara.z_list[i + 1]
        iterpara.A2[i] = iterpara.A2[i + 1] + q[j] * exp(- 1.0im * (k_x * x[j] + k_y * y[j]) - s * para.α * z[j])
    end

    return nothing
end

# notice that B have nothing to do with (s, w), so that it only needed to be updated once for each K
function update_iterpara_B12!(iterpara::IterPara, q::Vector{T}, x::Vector{T}, y::Vector{T}, z::Vector{T}, para::SoEwald2DLongInteraction{T}, K::Tuple{T, T, T}) where{T<:Number}
    k_x, k_y, k = K

    #update B
    iterpara.B1[1] = zero(ComplexF64)
    iterpara.B2[para.n_atoms] = zero(ComplexF64)
    for i in 2:para.n_atoms
        iterpara.B1[i] = iterpara.B1[i - 1]
        for m in iterpara.m_list[i - 1] : iterpara.m_list[i] - 1
            j = iterpara.z_list[m]
            iterpara.B1[i] += 2 * q[j] * exp( - 1.0im * (k_x * x[j] + k_y * y[j]) + k * z[j])
        end 
    end

    return nothing
end

function update_iterpara_C12!(iterpara::IterPara, q::Vector{T}, x::Vector{T}, y::Vector{T}, z::Vector{T}, para::SoEwald2DLongInteraction{T}, s::ComplexF64, K::Tuple{T, T, T}) where{T<:Number}
    k_x, k_y, k = K

    #update C
    iterpara.C[1] = zero(ComplexF64)
    for i in 2:para.n_atoms
        iterpara.C[i] = iterpara.C[i - 1]
        for m in iterpara.m_list[i - 1] : iterpara.m_list[i] - 1
            j = iterpara.z_list[m]
            iterpara.C[i] += q[j] * exp( - 1.0im * (k_x * x[j] + k_y * y[j]) + k * z[j] + s * para.α * z[j])
        end 
    end

    return nothing
end

function update_iterpara_D12!(iterpara::IterPara, q::Vector{T}, x::Vector{T}, y::Vector{T}, z::Vector{T}, para::SoEwald2DLongInteraction{T}, s::ComplexF64, K::Tuple{T, T, T}) where{T<:Number}
    k_x, k_y, k = K
    n_atoms = para.n_atoms
    z_list = iterpara.z_list
    m_list = iterpara.m_list

    # n is the max value for the summation
    # n = para.α * para.L[3] * s
    if real(k - s * para.α) < 0
        # use the reverse mode iter
        iterpara.D[n_atoms] = zero(ComplexF64)
        for i in n_atoms - 1:-1:m_list[n_atoms]
            j = z_list[i]
            iterpara.D[n_atoms] += q[j] * exp( - 1.0im * (k_x * x[j] + k_y * y[j]) + k * z[j] - s * para.α * z[j])
        end

        for i in n_atoms - 1:-1:1
            li = z_list[i]
            iterpara.D[i] = iterpara.D[i + 1] - q[li] * exp(-1.0im * (k_x * x[li] + k_y * y[li]) + k * z[li] - s * para.α * z[li])
            for j in m_list[i + 1] - 1 : -1 : m_list[i]
                lj = z_list[j]
                iterpara.D[i] += q[lj] * exp(-1.0im * (k_x * x[lj] + k_y * y[lj]) + k * z[lj] - s * para.α * z[lj])
            end
        end
    else
        # use the forward mode iter
        iterpara.D[1] = zero(ComplexF64)
        for i in 2:n_atoms
            iterpara.D[i] = iterpara.D[i - 1]
            for j in m_list[i - 1]:m_list[i] - 1
                lj = z_list[j]
                iterpara.D[i] -= q[lj] * exp(-1.0im * (k_x * x[lj] + k_y * y[lj]) + k * z[lj] - s * para.α * z[lj])
            end
            li = z_list[i - 1]
            iterpara.D[i] += q[li] * exp(-1.0im * (k_x * x[li] + k_y * y[li]) + k * z[li] - s * para.α * z[li])
        end
    end
    
    return nothing
end

function force_sum_k0!(q::Vector{T}, z::Vector{T}, para::SoEwald2DLongInteraction{T}, soepara::SoePara{ComplexF64}, iterpara::IterPara, sum::Vector{Point{3, ComplexF64}}) where{T}

    α = para.α
    z_list = iterpara.z_list
    n_atoms = para.n_atoms

    iterpara.A[1] = zero(ComplexF64)
    iterpara.B[n_atoms] = zero(ComplexF64)
    iterpara.C[1] = zero(ComplexF64)
    iterpara.D[n_atoms] = zero(ComplexF64)

    # irrelevant to (s, w)
    for i in n_atoms - 1:-1:1
        l = z_list[i + 1]
        iterpara.B[i] = iterpara.B[i + 1] + q[l]
    end
    for i in 2:n_atoms
        l = z_list[i - 1]
        iterpara.A[i] = iterpara.A[i - 1] + q[l]
    end
    for i in 1:n_atoms
        l = z_list[i]
        q_i = q[l]
        sum[l] -= Point(zero(ComplexF64), zero(ComplexF64), ComplexF64(q_i * (iterpara.A[i] - iterpara.B[i]) / (4 *  para.L[1] * para.L[2])))
    end

    # relevant to (s, w)
    for (s, w) in soepara.sw
        for i in 2:n_atoms
            l = z_list[i - 1]
            iterpara.C[i] = iterpara.C[i - 1] + q[l] * exp(s * α * z[l])
        end
        for i in n_atoms - 1:-1:1
            l = z_list[i + 1]
            iterpara.D[i] = iterpara.D[i + 1] + q[l] * exp( - s * α * z[l])
        end
        for i in 1:n_atoms
            l = z_list[i]
            q_i = q[l]
            z_i = z[l]
            sum[l] += Point(zero(ComplexF64), zero(ComplexF64), w * q_i * (iterpara.C[i] * exp( - s * α * z_i) - iterpara.D[i] * exp(s * α * z_i)) / (4 *  para.L[1] * para.L[2]))
        end
    end
    return nothing
end

function force_sum_k!(K::Tuple{T, T, T}, q::Vector{T}, x::Vector{T}, y::Vector{T}, z::Vector{T}, para::SoEwald2DLongInteraction{T}, soepara::SoePara{ComplexF64}, iterpara::IterPara, sum::Vector{Point{3, ComplexF64}}) where{T}
    k_x, k_y, k = K
    α = para.α
    update_iterpara_m!(iterpara, z, k / (2 * α^2))
    
    # compute the part irrelevant to SOE
    update_iterpara_B!(iterpara, q, x, y, z, para, K)

    for i in 1:para.n_atoms
        j = iterpara.z_list[i]
        sum_xy = q[j] * exp(1.0im * (k_x * x[j] + k_y * y[j]) - k * z[j]) * iterpara.B[i]
        sum[j] += Point(1.0im * k_x * sum_xy / k, 1.0im * k_y * sum_xy / k, zero(ComplexF64)) / (8 * para.L[1] * para.L[2])
    end

    #compute the SOE part
    for (s, w) in soepara.sw
        #update A, C and D
        update_iterpara_A12!(iterpara, q, x, y, z, para, s, K)
        update_iterpara_C!(iterpara, q, x, y, z, para, s, K)
        update_iterpara_D!(iterpara, q, x, y, z, para, s, K)
        
        for i in 1:para.n_atoms
            j = iterpara.z_list[i]
            sum_k_soe_1 = q[j] * exp(1.0im * (k_x * x[j] + k_y * y[j])) * (exp(- s*para.α*z[j]) * iterpara.A1[i] + exp(s*para.α*z[j]) * iterpara.A2[i])
            sum_k_soe_2 = - q[j] * exp(1.0im * (k_x * x[j] + k_y * y[j]) - k * z[j] - s*para.α*z[j]) * iterpara.C[i]
            sum_k_soe_3 = q[j] * exp(1.0im * (k_x * x[j] + k_y * y[j]) - k * z[j] + s*para.α*z[j]) * iterpara.D[i]

            sum_xy_k = s * w * exp(-k^2/(4 * para.α^2)) / (s * para.α + k) * para.α * sum_k_soe_1 + w * (exp(s * k / (2 * para.α)) * sum_k_soe_2 + exp(- s * k / (2 * para.α)) * sum_k_soe_3)

            sum[j] += Point(1.0im * k_x * sum_xy_k / k, 1.0im * k_y * sum_xy_k / k, zero(ComplexF64)) / (8 * para.L[1] * para.L[2])
        end
    end

    return nothing
end

function SoEwald2D_Fl!(interaction::SoEwald2DLongInteraction{T}, sys::MDSys, info::SimulationInfo{T}) where{T<:Number}

    revise_interaction!(interaction, sys, info)

    iterpara = interaction.iterpara
    soepara = interaction.soepara
    q = interaction.q
    x = interaction.x
    y = interaction.y
    z = interaction.z
    mass = interaction.mass

    update_iterpara_z!(iterpara, z)

    sum_temp = [Point(3, zero(ComplexF64)) for _=1:interaction.n_atoms]
    force_sum_k0!(q, z, interaction, soepara, iterpara, sum_temp)

    for i in 1:size(interaction.k_set, 1)
        K = interaction.k_set[i]
        force_sum_k!(K, q, x, y, z, interaction, soepara, iterpara, sum_temp)
    end

    info.acceleration .+= real.(sum_temp) ./ mass

    return nothing 
end