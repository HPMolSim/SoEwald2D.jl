# Nearest periodic image for a slab geometry, replacing this package's use of
# ExTinyMD's `position_check3D` with a `Q2dBoundary`.
#
# What `position_check3D` actually did here, since that is what has to be
# reproduced: `Q2dBoundary(Lx, Ly, Lz)` is `Boundary((Lx,Ly,Lz), (1,1,0))`, so
# its `mz` loop runs only at `mz = 0`. It therefore scanned `mx, my in -1:1`,
# wrapped x and y, left z alone, and returned the FIRST image whose **full 3-D**
# squared distance was below the cutoff -- plus an all-zero sentinel triple when
# none was, which forced every caller to guard with `iszero(r_sq)`.
#
# `_min_image_slab` is the true minimum image for any cutoff (no `-1:1` scan
# limit, so it is also right for a coordinate that has drifted outside the unit
# cell), and callers test the returned `r_sq` against `r_c^2` explicitly rather
# than against a sentinel.
#
# OPERATION ORDER IS LOAD-BEARING. The nearest in-plane image is
# `_wrap(dx, L) = dx - L*round(dx/L)` mathematically, but computing it that way
# and then reconstructing `coord_i` as `pos_j + dx` does NOT reproduce
# `position_check3D`'s arithmetic, which shifts `pos_i` by a whole number of
# periods and subtracts afterwards -- `(a - L) - b` and `(a - b) - L` differ in
# the last bits. Measured on three configurations, the `pos_j + dx` form moved
# the short-range energy by up to 2.2e-16 relative and the short-range forces by
# up to 1.5e-14 relative. So the image count is taken first and applied to
# `pos_i`, exactly as `coord_1 + Point(mx*Lx, my*Ly, 0)` did, which makes the
# whole replacement bit-identical.
@inline _image_shift(d::T, L::T) where {T} = -round(d / L) * L

"""
    _min_image_slab(pos_i, pos_j, L) -> (coord_i, coord_j, r_sq)

Slab-geometry nearest image: `x` and `y` wrap under `L[1]`/`L[2]`, `z` is a plain
difference (the slab axis is not periodic). Returns `pos_i` shifted to its nearest
in-plane image of `pos_j`, `pos_j` unchanged, and the **full three-dimensional**
squared distance `dx^2 + dy^2 + dz^2` between them.

The 3-D distance is the point of this helper and the reason it is not QuasiEwald's
`_min_image_q2d`, which returns the in-plane distance instead: SoEwald2D's real-space
sum is `erfc(α·r)/r` in the true separation `r`, and its cutoff test is on `r`, not on
the in-plane `ρ`. Swapping in an in-plane version would silently change every
short-range pair's cutoff test and its `erfc` argument -- the same class of mistake
that made `Ewald2D` + `CellListQ2D` return `+0.0238` against a true `−0.1539` in
Phase 1 of this project.

`pos_i`/`pos_j` need only support `p[1]`/`p[2]`/`p[3]` indexing, so an
`SVector{3,T}`, an `NTuple{3,T}` and ExTinyMD's `Point{3,T}` all work with no
conversion layer.

Only correct as the *whole* short-range story while `r_c < min(Lx, Ly) / 2`; at or
above half the box a second image is also inside the cutoff and the single nearest
one is not enough. [`SoEwald2DShortPlan`](@ref) enforces that bound at construction.
"""
@inline function _min_image_slab(pos_i, pos_j, L::NTuple{3, T}) where {T}
    x_i = T(pos_i[1]); y_i = T(pos_i[2]); z_i = T(pos_i[3])
    x_j = T(pos_j[1]); y_j = T(pos_j[2]); z_j = T(pos_j[3])

    coord_i = SVector{3, T}(x_i + _image_shift(x_i - x_j, L[1]),
                            y_i + _image_shift(y_i - y_j, L[2]),
                            z_i)
    coord_j = SVector{3, T}(x_j, y_j, z_j)

    dx = coord_i[1] - coord_j[1]
    dy = coord_i[2] - coord_j[2]
    dz = coord_i[3] - coord_j[3]
    r_sq = dx^2 + dy^2 + dz^2

    return coord_i, coord_j, r_sq
end

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

# ============================================================================
# Framework-free plans.
#
# `SoEwald2DShortInteraction`/`SoEwald2DLongInteraction` used to be
# `ExTinyMD.AbstractInteraction` structs holding both the method's parameters
# and the MD side's gathered state. The plans below are the framework-free
# replacement: the same physics, constructed and queried from plain arrays via
# `SoEwald2D.energy`/`force`/`force!` (defined alongside the rest of the
# energy/force machinery in energy/energy_short.jl, energy/energy_long.jl,
# force/force_short.jl, force/force_long.jl).
#
# `mass` and `acceleration` are gone from the long plan. A framework-free
# solver returns a FORCE and leaves mass-division to the caller, exactly as
# ExTinyMD's own electrostatics adapter does; the ExTinyMD wrapper's
# `update_acceleration!` does that division. Verified bitwise, not assumed:
# the old `SoEwald2D_Fl!` did `acceleration -= Point(Fx[i], Fy[i], Fz[i]) /
# mass[i]`, ExTinyMD's `Point` `/` and `-` are both plain componentwise
# operations (`Point(y.coo ./ x)`, `Point(x.coo .- y.coo)`), and IEEE-754 makes
# `a - (f/m)` and `a + ((-f)/m)` the same double for every input -- negation is
# exact and division is correctly rounded, so `(-f)/m == -(f/m)` bit for bit.
# The force buffer therefore stores `-Fx[i]` and the wrapper adds `f/m`, with
# no change to any computed value.
# ============================================================================

"""
    SoEwald2DShortPlan(ϵ_0, L, s, α, n_atoms, r_c)

Framework-free short-range (real-space) plan for the SOE Ewald2D method: pure
parameters, no ExTinyMD dependency, nothing MD-specific. Query with
[`SoEwald2D.energy`](@ref), [`SoEwald2D.force`](@ref) or
[`SoEwald2D.force!`](@ref) against plain array-of-structs positions
(`Vector{SVector{3,T}}` canonical, but anything supporting `p[1]`/`p[2]`/`p[3]`
indexing works, including `NTuple{3,T}` and ExTinyMD's `Point{3,T}`) and a
plain charge vector.

`r_c` must satisfy `r_c < min(Lx, Ly) / 2`; anything else throws an
`ArgumentError` (the short-range sum uses the single nearest in-plane periodic
image, which is only the whole story below half the box).
"""
struct SoEwald2DShortPlan{T}
    ϵ_0::T
    L::NTuple{3, T}
    s::T
    α::T
    n_atoms::Int64

    r_c::T
end

function SoEwald2DShortPlan(ϵ_0::T, L::NTuple{3, T}, s::T, α::T, n_atoms::Int64, r_c::T) where{T<:Number}
    # `_min_image_slab` returns the single nearest in-plane image, which is the
    # only image inside the cutoff exactly when `r_c < min(Lx, Ly) / 2`. At or
    # beyond half the box a second image is also within `r_c` and the
    # short-range sum silently omits it: no exception, no warning, just a wrong
    # number. This is the only place that can catch it for a standalone
    # caller, since nothing else in this package ever looks at the unit cell.
    if !(r_c < min(L[1], L[2]) / 2)
        throw(ArgumentError(
            "SoEwald2DShortPlan requires r_c < min(Lx, Ly) / 2, but got " *
            "r_c = $r_c with (Lx, Ly) = ($(L[1]), $(L[2])), i.e. " *
            "min(Lx, Ly) / 2 = $(min(L[1], L[2]) / 2). The real-space sum uses " *
            "the single nearest in-plane periodic image, which is only the whole " *
            "story below half the box; at or above it a second image is also " *
            "inside the cutoff and is silently dropped. Reduce r_c (equivalently " *
            "reduce s = α·r_c), or enlarge Lx/Ly."))
    end
    return SoEwald2DShortPlan{T}(ϵ_0, L, s, α, n_atoms, r_c)
end

"""
    SoEwald2DLongPlan(ϵ_0, L, s, α, n_atoms, k_c, soepara; rbm = false, rbm_p = 0, parallel = true, rng = MersenneTwister(123))

Framework-free long-range (reciprocal-space) plan for the SOE Ewald2D method.
Same parameters as the old `SoEwald2DLongInteraction`, minus `mass` and
`acceleration`: [`SoEwald2D.force!`](@ref) returns a force, not an
acceleration, so this plan carries no notion of mass at all.

The kernels below `energy_sum!`/`force_sum` are structure-of-arrays, so the
plan owns `q`, `x`, `y`, `z` as scratch and every query scatters the caller's
array-of-structs `poses`/`charges` into them (see `_scatter_long!`). The
caller's arrays are never written to.

Query with [`SoEwald2D.energy`](@ref)/[`SoEwald2D.force`](@ref)/
[`SoEwald2D.force!`](@ref).
"""
struct SoEwald2DLongPlan{T}
    ϵ_0::T
    L::NTuple{3, T}
    s::T
    α::T
    n_atoms::Int64

    k_c::T
    k_set::Vector{Tuple{T, T, T}}
    soepara::SoePara{ComplexF64}
    rbm::Bool # whether to use random batch method
    rbm_p::Int
    P::T
    prob::ProbabilityWeights{T}
    indice::Vector{Int}

    # Plan-owned structure-of-arrays scratch, refilled from the caller's AoS
    # arrays by `_scatter_long!` at the start of every query.
    q::Vector{T}
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
    iterpara::IterPara
    adpara::AdPara

    parallel::Bool
    rng
end

function SoEwald2DLongPlan(ϵ_0::T, L::NTuple{3, T}, s::T, α::T, n_atoms::Int64, k_c::T, soepara::SoePara{ComplexF64}; rbm::Bool = false, rbm_p::Int=0, parallel::Bool = true, rng = MersenneTwister(123)) where{T<:Number}

    k_set = Vector{Tuple{T, T, T}}()
    if rbm == false
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
        P = zero(T)
        prob = ProbabilityWeights([zero(T)])
    else
        k_set, P, prob = generate_K_set(α, L, k_c)
    end

    indice = zeros(Int, rbm_p)

    q = zeros(T, n_atoms)
    x = zeros(T, n_atoms)
    y = zeros(T, n_atoms)
    z = zeros(T, n_atoms)

    iterpara = IterPara(n_atoms)
    adpara = AdPara(n_atoms)

    return SoEwald2DLongPlan(ϵ_0, L, s, α, n_atoms, k_c, k_set, soepara, rbm, rbm_p, P, prob, indice, q, x, y, z, iterpara, adpara, parallel, rng)
end

"""
    _scatter_long!(plan, poses, charges) -> nothing

Fill the long plan's structure-of-arrays scratch from array-of-structs
`poses`/`charges`. This is the framework-free replacement for the old
`revise_interaction!(interaction, sys, ExTinyMD.SimulationInfo)`, which read
the same values out of `sys.atoms`/`info.particle_info` instead. Neither
argument is mutated -- `poses` and `charges` are read only, and everything
written lives on the plan.

`poses[i]` need only support `p[1]`/`p[2]`/`p[3]` indexing.
"""
function _scatter_long!(plan::SoEwald2DLongPlan{T}, poses, charges) where{T}
    @inbounds for i in 1:plan.n_atoms
        plan.q[i] = charges[i]
        p = poses[i]
        plan.x[i] = p[1]
        plan.y[i] = p[2]
        plan.z[i] = p[3]
    end
    return nothing
end
