function energy_sum_k_S1(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, Int64}, soepara::SoePara{ComplexF64}, iterpara::IterPara) where{T<:Number}
    k_x, k_y, k = K
    α = para.α
    update_iterpara_z!(iterpara, z)
    update_iterpara_m!(iterpara, z, k / (2 * α^2))
    
    sum_k = zero(ComplexF64)

    #compute the SOE part
    for (s, w) in soepara.sw
        #update A, C and D
        sum_k_soe = zero(ComplexF64)
        update_iterpara_A!(iterpara, q, x, y, z, para, s, K)
        for i in 1:para.n_atoms
            j = iterpara.z_list[i]
            sum_k_soe += q[j] * exp(1.0im * (k_x * x[j] + k_y * y[j]) + k * z[j] - s*para.α*z[j]) * iterpara.A[i]
        end

        sum_k += w * sum_k_soe * exp(- s * k / (2 * para.α))
    end

    return 2/k * (sum_k)
end

function energy_sum_k_S2(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, Int64}, soepara::SoePara{ComplexF64}, iterpara::IterPara) where{T<:Number}
    k_x, k_y, k = K
    α = para.α
    update_iterpara_z!(iterpara, z)
    update_iterpara_m!(iterpara, z, k / (2 * α^2))

    # compute the part irrelevant to SOE
    update_iterpara_B!(iterpara, q, x, y, z, para, K)
    sum_k0 = zero(ComplexF64)
    for i in 1:para.n_atoms
        # sum_k += q[i] * q[i] * erfc(k/(2 * para.α))
        j = iterpara.z_list[i]
        sum_k0 += q[j] * exp(1.0im * (k_x * x[j] + k_y * y[j]) - k * z[j]) * iterpara.B[i]
    end
    

    # compute the SOE part
    sum_k1 = zero(ComplexF64)
    sum_k2 = zero(ComplexF64)
    for (s, w) in soepara.sw
        #update A, C and D
        sum_k_soe1 = zero(ComplexF64)
        sum_k_soe2 = zero(ComplexF64)
        update_iterpara_C!(iterpara, q, x, y, z, para, s, K)
        update_iterpara_D!(iterpara, q, x, y, z, para, s, K)

        for i in 1:para.n_atoms
            j = iterpara.z_list[i]
            sum_k_soe1 += - q[j] * exp(1.0im * (k_x * x[j] + k_y * y[j]) - k * z[j] - s*para.α*z[j]) * iterpara.C[i]
            sum_k_soe2 += q[j] * (exp(1.0im * (k_x * x[j] + k_y * y[j]) - k * z[j] + s*para.α*z[j]) * iterpara.D[i])
        end

        sum_k1 += w * exp(s * k / (2 * para.α)) * sum_k_soe1
        sum_k2 += w * exp(- s * k / (2 * para.α)) * sum_k_soe2
    end
    return 2/k * (sum_k0), 2/k * (sum_k1), 2/k * (sum_k2)
end

function big_iterpara_A(q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, Int64}, s::ComplexF64, K::Tuple{T, T, T}) where{T<:Number}
    k_x, k_y, k = K

    A = [big(zero(ComplexF64)) for _=1:para.n_atoms]
    z_list = sortperm(z)
    for i in 2:para.n_atoms
        for l in 1:i-1
            j = z_list[l]
            A[i] += q[j] * exp(big(- 1.0im * (k_x * x[j] + k_y * y[j]) - k * z[j] + s * para.α * z[j]))
        end
    end

    return A
end

function big_sum_A(q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, para::SoEwald2DPara{T, Int64}, soepara::SoePara{ComplexF64}, K::Tuple{T, T, T}) where{T<:Number}
    k_x, k_y, k = K

    result_1 = big(zero(ComplexF64))
    result_2 = big(zero(ComplexF64))
    result_3 = big(zero(ComplexF64))
    result_4 = big(zero(ComplexF64))

    n = k * para.L[3]
    z_list = sortperm(z)
    for (s, w) in soepara.sw
        for i in 1:para.n_atoms
            A = big(zero(ComplexF64))
            for l in 1:i-1
                j = z_list[l]
                A += q[j] * exp(big(- 1.0im * (k_x * x[j] + k_y * y[j]) - k * z[j] + s * para.α * z[j]))
            end
            j = z_list[i]
            result_1 += q[j] * exp(big(log(w) + 1.0im * (k_x * x[j] + k_y * y[j]) + k * z[j] - s * para.α * z[j] - s * k / (2 * para.α))) * A

            for j in 1:i - 1
                result_2 += q[i] * q[j] * exp(big(log(w) + 1.0im * (k_x * (x[i] - x[j]) + k_y * (y[i] - y[j])) + k * abs(z[i] - z[j]) - s * para.α * abs(z[i] - z[j]) - s * k / (2 * para.α)))
            end
        end
    end

    e_1 = []
    e_2 = []
    t = []
    for i in 1:para.n_atoms
        for j in 1:i - 1
            t1 = q[i] * q[j] * exp(big( 1.0im * (k_x * (x[i] - x[j]) + k_y * (y[i] - y[j])) + k * abs(z[i] - z[j])))
            result_4 += t1 * erfc(big(k/(2*para.α) + para.α * abs(z[i] - z[j])))
            temp = zero(ComplexF64)

            s, w = soepara.sw[1]
            temp = (w)^(1/n) * exp(big(- s * para.α * abs(z[i] - z[j]) - s * k / (2 * para.α)) / n)
            for l in 2:size(soepara.sw)[1]
                s, w = s = soepara.sw[l]
                a = exp(big(- s * para.α * abs(z[i] - z[j]) - s * k / (2 * para.α)) / n)
                if abs(a) > abs(temp)
                    temp = a * (w + (temp / a)^n)^(1/n)
                else
                    temp = temp * (w * (a / temp)^n + 1)^(1/n)
                end
            end

            t2 =  exp(big( 1.0im * (k_x * (x[i] - x[j]) + k_y * (y[i] - y[j])) + k * abs(z[i] - z[j])) / n)
            result_3 += q[i] * q[j] * (t2 * temp)^n
            if abs(t1 * temp - t1 * erfc(big(k/(2*para.α) + para.α * abs(z[i] - z[j])))) > 1e-5
                @show t1, temp - erfc(big(k/(2*para.α) + para.α * abs(z[i] - z[j])))
            end

            a = temp - erfc(big(k/(2*para.α) + para.α * abs(z[i] - z[j])))
            b = result_3 - result_4
            push!(e_1, a)
            push!(e_2, b)
            push!(t, t1 * temp - t1 * erfc(big(k/(2*para.α) + para.α * abs(z[i] - z[j]))) )
        end
    end

    return 2/k * result_1, 2/k * result_2, 2/k * result_3, 2/k * result_4, e_1, e_2, t
end