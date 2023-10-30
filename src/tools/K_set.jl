function generate_K_set(α::T, L::NTuple{3, T}, set_size::Int) where{T <: Number}

    L_x, L_y, L_z = L
    k_c = sqrt(64 * α^2 * log(10)) + 4π / sqrt(L_x * L_y)

    mx_max = ceil(Int, k_c * L_x / 2π) + 1
    my_max = ceil(Int, k_c * L_y / 2π) + 1

    K_set = Vector{NTuple{3, T}}()
    Prob = Vector{T}()

    sum_K = zero(T)

    for m_x in - mx_max : mx_max
        for m_y in - my_max : my_max
            k_x = m_x * 2π / L_x
            k_y = m_y * 2π / L_y
            k = sqrt(k_x^2 + k_y^2)
            if k != 0 && k < k_c
                push!(K_set, (k_x, k_y, k))
                prob = exp(-k^2 / (4 * α^2))
                push!(Prob, prob)
                sum_K += prob
            end
        end
    end
    Prob ./= sum_K
    K_set = sample(K_set, ProbabilityWeights(Prob), set_size)

    return K_set, sum_K
end

"""
function vec_divider(a::Vector{T}, n::Int64) 
Divide an array into n parts with almost equal length, if n > length(a);
else it will result length(a) parts.
Used for parallel computing.

"""
function vec_divider(a::Vector{T}, n::Int64) where{T}
    if n == 0
        # @warn "Input n is 0, automatically set to 1."
        n = 1
    end
    da = [Vector{T}() for i in 1:min(n, length(a))]
    for i in 1:length(a)
        id = i % n
        if id == 0
            id = n
        end
        push!(da[id], a[i])
    end
    return da
end