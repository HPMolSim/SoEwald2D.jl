function generate_K_set(α::T, L::NTuple{3, T}, k_c::T, set_size::Int) where{T <: Number}

    L_x, L_y, L_z = L

    mx_max = ceil(Int, k_c * L_x / 2π) + 1
    my_max = ceil(Int, k_c * L_y / 2π) + 1

    K_set = Vector{NTuple{3, T}}()
    Prob = Vector{T}()

    sum_K = big(0.0)

    for m_x in - mx_max : mx_max
        for m_y in - my_max : my_max
            k_x = m_x * 2π / L_x
            k_y = m_y * 2π / L_y
            k = sqrt(k_x^2 + k_y^2)
            if 0 < k <= k_c
                push!(K_set, (k_x, k_y, k))
                prob = exp(-k^2 / (4 * α^2))
                push!(Prob, prob)
                sum_K += big(prob)
            end
        end
    end
    Prob ./= T(sum_K)
    K_set = sample(K_set, ProbabilityWeights(Prob), set_size)

    return K_set, T(sum_K)
end