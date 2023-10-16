struct SoePara{T} 
    sw::Vector{Tuple{T, T}}
end

mutable struct IterPara
    A::Vector{ComplexF64}
    B::Vector{ComplexF64}
    z_list::Vector{Int64}
end

# due to the soe approach is defined in complex space, we directlly defined the vector as ComplexF64
function IterPara(n_atoms::Int64)
    A = zeros(ComplexF64, n_atoms)
    B = zeros(ComplexF64, n_atoms)

    z_list = zeros(Int64, n_atoms)

    return IterPara(A, B, z_list)
end

# this function is possibly not needed 
function revise_iterpara!(iterpara::IterPara)
    n_atoms = length(iterpara.A)
    for i = 1:n_atoms
        iterpara.A[i] = zero(ComplexF64)
        iterpara.B[i] = zero(ComplexF64)
        iterpara.z_list[i] = zero(Int64)
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

struct SoEwald2DLongInteraction{T} <: ExTinyMD.AbstractInteraction
    ϵ_0::T
    L::NTuple{3, T}
    accuracy::T
    α::T
    n_atoms::Int64

    k_c::T
    k_set::Vector{Tuple{T, T, T}}
    soepara::SoePara{ComplexF64}

    q::Vector{T}
    mass::Vector{T}
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
    acceleration::Vector{Point{3, T}}
    iterpara::IterPara
    adpara::AdPara
end

function SoEwald2DLongInteraction(ϵ_0::T, L::NTuple{3, T}, accuracy::T, α::T, n_atoms::Int64, k_c::T, soepara::SoePara{ComplexF64}) where{T<:Number}

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

    q = zeros(T, n_atoms)
    mass = zeros(T, n_atoms)
    x = zeros(T, n_atoms)
    y = zeros(T, n_atoms)
    z = zeros(T, n_atoms)
    acceleration = Vector{Point{3, T}}(undef, n_atoms)

    iterpara = IterPara(n_atoms)
    adpara = AdPara(n_atoms)

    return SoEwald2DLongInteraction(ϵ_0, L, accuracy, α, n_atoms, k_c, k_set, soepara, q, mass, x, y, z, acceleration, iterpara, adpara)
end

struct SoEwald2DShortInteraction{T} <: ExTinyMD.AbstractInteraction
    ϵ_0::T
    L::NTuple{3, T}
    accuracy::T
    α::T
    n_atoms::Int64

    r_c::T
end

function SoEwald2DShortInteraction(ϵ_0::T, L::NTuple{3, T}, accuracy::T, α::T, n_atoms::Int64, r_c::T) where{T<:Number}
    return SoEwald2DShortInteraction{T}(ϵ_0, L, accuracy, α, n_atoms, r_c)
end

function revise_interaction!(interaction::SoEwald2DLongInteraction{T}, sys::MDSys{T}, info::SimulationInfo{T}) where{T<:Number}

    for i in 1:length(interaction.q)
        interaction.q[i] = sys.atoms[info.particle_info[i].id].charge
        interaction.mass[i] = sys.atoms[info.particle_info[i].id].mass
        interaction.x[i], interaction.y[i], interaction.z[i] = info.particle_info[i].position
    end

    return nothing
end