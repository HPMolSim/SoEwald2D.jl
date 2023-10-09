export SoePara, SoEwald2DPara, AdPara, IterPara, revise_adpara!, update_iterpara!, SoEwald2DLongInteraction, revise_interaction!

struct SoePara{T} 
    sw::Vector{Tuple{T, T}}
end

# the default 18 terms approach
SoePara() = SoePara{ComplexF64}([
    (4.56160429581067 - 8.10664614462486im, 1.7841998933545696e-7 + 7.343003636293712e-7im),
    (4.56160429581067 + 8.10664614462486im, 1.7841998935002611e-7 - 7.343003636293118e-7im),
    (4.67614751620509 - 6.66817205403247im, -9.556499453539262e-5 - 0.00015001811262319774im),
    (4.67614751620509 + 6.66817205403247im, -9.556499453567014e-5 + 0.00015001811262314695im),
    (4.75706771210077 - 5.46277771275372im, 0.007705869731927555 + 0.002956305642083839im),
    (4.75706771210077 + 5.46277771275372im, 0.007705869731931252 - 0.0029563056420808262im),
    (4.81739468471801 - 4.37139574323745im, -0.13131561540899633 + 0.0581715387810661im),
    (4.81739468471801 + 4.37139574323745im, -0.13131561540901635 - 0.0581715387810941im),
    (4.86201933938979 - 3.34750236244862im, 0.3353050941709935 - 1.1598168978326757im),
    (4.86201933938979 + 3.34750236244862im, 0.3353050941710323 + 1.15981689783286im),
    (4.89345823682377 - 2.3658861273834im, 3.3749677449815447 + 4.418684212942269im),
    (4.89345823682377 + 2.3658861273834im, 3.3749677449817934 - 4.41868421294285im),
    (4.91338156910332 - 1.41021110952036im, -14.965174561023547 + 0.8807329873608107im),
    (4.91338156910332 + 1.41021110952036im, -14.965174561024503 - 0.8807329873602142im),
    (4.92300187934591 - 0.468589313154398im, 11.878606854122312 - 21.3765363724229im),
    (4.92300187934591 + 0.468589313154398im, 11.878606854123603 + 21.376536372422756im),
    (0.0334524767879181 - 0.0245907245253121im, 4.394108031340191e-15 + 6.7735884437017495e-15im),
    (0.0334524767879181 + 0.0245907245253121im, 4.394104404590236e-15 - 6.773593136580561e-15im)
])

struct SoExPara{T} 
    sw::Vector{Tuple{T, T}}
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

    return SoEwald2DLongInteraction(ϵ_0, L, accuracy, α, n_atoms, k_c, k_set, soepara, q, mass, x, y, z, acceleration, iterpara)
end

struct SoEwald2DShortInteraction{T}
    ϵ_0::T
    L::NTuple{3, T}
    accuracy::T
    α::T
    n_atoms::Int64
end

function revise_interaction!(interaction::SoEwald2DLongInteraction{T}, sys::MDSys{T}, info::SimulationInfo{T}) where{T<:Number}

    for i in 1:length(interaction.q)
        interaction.q[i] = sys.atoms[info.particle_info[i].id].charge
        interaction.mass[i] = sys.atoms[info.particle_info[i].id].mass
        interaction.x[i], interaction.y[i], interaction.z[i] = info.particle_info[i].position
    end

    return nothing
end