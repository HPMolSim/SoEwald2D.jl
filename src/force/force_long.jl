function energy_sum_k!(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, n_atoms::Int64, α::T, soepara::SoePara{ComplexF64}, iterpara::IterPara, U::Array{T}) where{T<:Number}
    U[1] = energy_sum_k(K, q, x, y, z, n_atoms, α, soepara, iterpara)
    return nothing
end

function energy_sum_k0!(q::Array{T}, z::Array{T}, n_atoms::Int64, α::T, soepara::SoePara{ComplexF64}, iterpara::IterPara, U::Array{T}) where{T<:Number}
    U[1] = energy_sum_k0(q, z, n_atoms, α, soepara, iterpara)
    return nothing
end

function force_sum_k(K::Tuple{T, T, T}, q::Array{T}, x::Array{T}, y::Array{T}, z::Array{T}, n_atoms::Int64, α::T, soepara::SoePara{ComplexF64}, iterpara::IterPara, adpara::AdPara) where{T<:Number}

    revise_adpara!(adpara, n_atoms)

    autodiff(ReverseWithPrimal, energy_sum_k!, Const(K), Const(q), Duplicated(x, adpara.Fx), Duplicated(y, adpara.Fy), Duplicated(z, adpara.Fz), Const(n_atoms), Const(α), Const(soepara), Duplicated(iterpara, adpara.iterpara_t), Duplicated(adpara.U, adpara.dU))

    return [adpara.Fx, adpara.Fy, adpara.Fz]
end

function force_sum_k0(q::Array{T}, z::Array{T}, n_atoms::Int64, α::T, soepara::SoePara{ComplexF64}, iterpara::IterPara, adpara::AdPara) where{T<:Number}

    revise_adpara!(adpara, n_atoms)
    autodiff(ReverseWithPrimal, energy_sum_k0!, Const(q), Duplicated(z, adpara.Fz), Const(n_atoms), Const(α), Const(soepara), Duplicated(iterpara, adpara.iterpara_t), Duplicated(adpara.U, adpara.dU))

    return [adpara.Fx, adpara.Fy, adpara.Fz]
end

# Returns the energy GRADIENT, not the force -- `[dU/dx, dU/dy, dU/dz]`, one
# array per axis. `force!` below negates it. This is the sign the
# pre-decoupling code carried too: `SoEwald2D_Fl!` did
# `acceleration -= Point(Fx, Fy, Fz) / mass`.
function force_sum(plan::SoEwald2DLongPlan{T}) where{T}

    iterpara = plan.iterpara
    soepara = plan.soepara
    adpara = plan.adpara
    q = plan.q
    x = plan.x
    y = plan.y
    z = plan.z

    rbm = plan.rbm
    rbm_p = plan.rbm_p
    P = plan.P
    prob = plan.prob

    α = plan.α
    n_atoms = plan.n_atoms
    L = plan.L
    k_set = plan.k_set
    ϵ_0 = plan.ϵ_0
    parallel = plan.parallel
    rng = plan.rng
    
    update_iterpara_z!(iterpara, z)

    F_k0 = - force_sum_k0(q, z, n_atoms, α, soepara, iterpara, adpara) * π / (L[1] * L[2])
    F_k = [zeros(T, n_atoms), zeros(T, n_atoms), zeros(T, n_atoms)]

    if parallel
        if rbm == false
            F_k = @distributed (+) for k in k_set
                exp(- k[3]^2 / (4 * α^2)) * force_sum_k(k, q, x, y, z, n_atoms, α, soepara, iterpara, adpara)
            end
        else
            random_k = sample(rng, k_set, prob, rbm_p)
            F_k = @distributed (+) for i in 1:rbm_p
                P / rbm_p * force_sum_k(random_k[i], q, x, y, z, n_atoms, α, soepara, iterpara, adpara)
            end
        end
    else
        if rbm == false
            for k in k_set
                F_k += exp(- k[3]^2 / (4 * α^2)) * force_sum_k(k, q, x, y, z, n_atoms, α, soepara, iterpara, adpara)
            end
        else
            random_k = sample(rng, k_set, prob, rbm_p)
            # `for i in indice` before the decoupling: `indice` is a *field* of
            # the plan and was never bound as a local here, so this branch --
            # rbm = true with parallel = false -- raised an unconditional
            # UndefVarError. `1:rbm_p` is what the parallel branch three lines
            # up uses and what `random_k` is indexed by.
            for i in 1:rbm_p
                F_k += P / rbm_p * force_sum_k(random_k[i], q, x, y, z, n_atoms, α, soepara, iterpara, adpara)
            end
        end
    end

    F_k *= π / (2 * L[1] * L[2])
    F_k .+= F_k0

    return F_k / (4π * ϵ_0)
end

"""
    SoEwald2D.force!(F, plan::SoEwald2DLongPlan, poses, charges) -> F
    SoEwald2D.force(plan::SoEwald2DLongPlan, poses, charges) -> Vector{SVector{3,T}}

Long-range (reciprocal-space) **force** from plain array-of-structs positions
and charges, written into `F` (filled, not accumulated into). Neither `poses`
nor `charges` is mutated -- the input is scattered into the plan's own
structure-of-arrays scratch first.

SIGN CONVENTION, AND IT IS A FLIP. [`force_sum`](@ref) returns the energy
*gradient*, and the pre-decoupling `SoEwald2D_Fl!` consumed it as one:
`acceleration -= Point(Fx[i], Fy[i], Fz[i]) / mass[i]`. This function returns
the FORCE, `-∇U`, so that the ExTinyMD wrapper accumulates with `+=` like every
other adapter in this family, and so that `force` means the same thing for both
plans. A finite difference of [`SoEwald2D.energy`](@ref) is what pins it down;
see test/plan.jl, which asserts the sign of every component rather than only
its magnitude.

No mass division is applied: a framework-free solver returns a force and leaves
that to the caller.

Not exported; call as `SoEwald2D.force!(...)`.
"""
function force!(F, plan::SoEwald2DLongPlan{T}, poses, charges) where{T<:Number}

    _scatter_long!(plan, poses, charges)
    revise_adpara!(plan.adpara, plan.n_atoms)

    Fx, Fy, Fz = force_sum(plan)

    @inbounds for i in 1:plan.n_atoms
        F[i] = SVector{3, T}(-Fx[i], -Fy[i], -Fz[i])
    end

    return F
end

"Allocating form of [`SoEwald2D.force!`](@ref)."
function force(plan::SoEwald2DLongPlan{T}, poses, charges) where{T<:Number}
    F = [SVector{3, T}(zero(T), zero(T), zero(T)) for _ in 1:plan.n_atoms]
    return force!(F, plan, poses, charges)
end
