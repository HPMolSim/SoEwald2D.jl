export SoEwald2DPara, AdPara, IterPara, revise_adpara!, update_iterpara!

# this structure contain all information needed for the calculation
struct SoEwald2DPara{T, TI}
    L::NTuple{3, T}
    r_c::T
    k_c::T
    α::T
    n_atoms::TI
    k_set::Vector{Tuple{T, T, T}}
end

function SoEwald2DPara(L::NTuple{3, T}, s::T, α::T, n_atoms::TI) where{T, TI}
    r_c = s * sqrt(α)
    k_c = 2*s/sqrt(α)
    k_set = Vector{Tuple{T, T, T}}()
    m_x = ceil(L[1] * k_c / 2π)
    m_y = ceil(L[2] * k_c / 2π)
    for i in - m_x : m_x
        for j in - m_y : m_y
            k_x = 2π * i / L[1]
            k_y = 2π * j / L[2]
            k = sqrt(k_x^2 + k_y^2)
            if k != 0 && k < k_c
                push!(k_set, (k_x, k_y, k))
            end
        end
    end
    return SoEwald2DPara{T, TI}(L, r_c, k_c, α, n_atoms, k_set)
end

mutable struct IterPara
    A::Vector{ComplexF64}
    B::Vector{ComplexF64}
    C::Vector{ComplexF64}
    D::Vector{ComplexF64}
    z_list::Vector{Int64}
    m_list::Vector{Int64}
end

# due to the soe approach is defined in complex space, we directlly defined the vector as ComplexF64
function IterPara(n_atoms::Int64)
    A = zeros(ComplexF64, n_atoms)
    B = zeros(ComplexF64, n_atoms)
    C = zeros(ComplexF64, n_atoms)
    D = zeros(ComplexF64, n_atoms)

    z_list = zeros(Int64, n_atoms)
    m_list = zeros(Int64, n_atoms)

    return IterPara(A, B, C, D, z_list, m_list)
end

# this function is possibly not needed 
function revise_iterpara!(iterpara::IterPara)
    n_atoms = length(iterpara.A)
    for i = 1:n_atoms
        iterpara.A[i] = zero(ComplexF64)
        iterpara.B[i] = zero(ComplexF64)
        iterpara.C[i] = zero(ComplexF64)
        iterpara.D[i] = zero(ComplexF64)
        iterpara.z_list[i] = zero(Int64)
        iterpara.m_list[i] = zero(Int64)
    end
    return nothing
end


mutable struct AdPara{T}
    Fx::Vector{T}
    Fy::Vector{T}
    Fz::Vector{T}
    U::Vector{T}
    dU::Vector{T}
    iterpara_t::IterPara
end

function AdPara(n_atoms::TI) where{TI<:Integer} 
    T = Float64
    Fx = zeros(T, n_atoms)
    Fy = zeros(T, n_atoms)
    Fz = zeros(T, n_atoms)
    U = [zero(T)]
    dU = [one(T)]
    iterpara_t = IterPara(n_atoms)
    return AdPara{T}(Fx, Fy, Fz, U, dU, iterpara_t)
end

function AdPara(T::DataType, n_atoms::TI) where{TI<:Integer} 
    Fx = zeros(T, n_atoms)
    Fy = zeros(T, n_atoms)
    Fz = zeros(T, n_atoms)
    U = [zero(T)]
    dU = [one(T)]
    iterpara_t = IterPara(n_atoms)
    return AdPara{T}(Fx, Fy, Fz, U, dU, iterpara_t)
end

function revise_adpara!(adpara::AdPara{T}, n_atoms::TI) where{T<:Number, TI<:Integer}
    for i in 1:n_atoms
        adpara.Fx[i] = zero(T)
        adpara.Fy[i] = zero(T)
        adpara.Fz[i] = zero(T)
    end
    adpara.U[1] = zero(T)
    adpara.dU[1] = one(T)
    adpara.iterpara_t = IterPara(n_atoms)
    return nothing
end