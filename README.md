# SoEwald2D

[![Build Status](https://github.com/ArrogantGao/SoEwald2D.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ArrogantGao/SoEwald2D.jl/actions/workflows/CI.yml?query=branch%3Amain)

`SoEwald2D.jl` is a `Julia` package for rapid calculation of the electrostatic
interaction in quasi-2D charged systems. It combines the SOE (sum-of-exponentials)
method with the Ewald2D method, so that its pairwise summation in k-space costs
$O(N)$ rather than the $O(N^2)$ of the original.

The package has two layers:

* a **framework-free core** — plan objects queried from plain array-of-structs
  positions and charges, with no MD framework anywhere; and
* a **thin MD adapter** in a package extension, so the same plans can drive
  `ExTinyMD.jl`'s `simulate!`.

`ExTinyMD` is a **weak** dependency: it is loaded only if you load it yourself.

## Getting started

```julia
pkg> add SoEwald2D               # core only
pkg> add SoEwald2D ExTinyMD      # core + the MD adapter
```

## Standalone usage (no ExTinyMD)

Build a plan, hand it positions and charges, get an energy or a force back.

```julia
using SoEwald2D, StaticArrays

n_atoms = 100
L = (20.0, 20.0, 20.0)
ϵ_0 = 1.0
α = 0.3
s = 2.0
r_c = s / α            # 6.667 < min(Lx, Ly) / 2 = 10.0
k_c = 2 * s * α

poses = [SVector(rand() * L[1], rand() * L[2], 0.5 + rand() * (L[3] - 1.0)) for _ in 1:n_atoms]
charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n_atoms]

short_plan = SoEwald2DShortPlan(ϵ_0, L, s, α, n_atoms, r_c)
long_plan  = SoEwald2DLongPlan(ϵ_0, L, s, α, n_atoms, k_c, SoePara())

E = SoEwald2D.energy(short_plan, poses, charges) +
    SoEwald2D.energy(long_plan,  poses, charges)

Fs = SoEwald2D.force(short_plan, poses, charges)     # allocates
Fl = SoEwald2D.force(long_plan,  poses, charges)

F = [SVector(0.0, 0.0, 0.0) for _ in 1:n_atoms]      # or fill a buffer
SoEwald2D.force!(F, long_plan, poses, charges)
```

Things worth knowing about that API:

* **`energy`, `force` and `force!` are deliberately not exported.** Call them
  qualified, as above. Several packages in this family define an `energy`, and
  exporting them all would make the name ambiguous and error on first use under
  `using EwaldSummations, SoEwald2D`. The plan *constructors* are exported,
  since those are what you need to discover.
* **`poses` is array-of-structs.** `Vector{SVector{3,T}}` is canonical, but the
  kernels index only `p[1]`, `p[2]` and `p[3]`, so a `Vector{NTuple{3,T}}` or a
  `Vector{ExTinyMD.Point{3,T}}` works identically, with no conversion layer and
  no numerical difference.
* **`r_c` must be strictly less than `min(Lx, Ly) / 2`.** `SoEwald2DShortPlan`
  throws an `ArgumentError` otherwise. The real-space sum uses the single
  nearest in-plane periodic image, which is the whole story only below half the
  box; at or above it a second image is also inside the cutoff and would be
  silently dropped.
* **Queries never mutate `poses` or `charges`.** The long plan scatters them
  into its own structure-of-arrays scratch, which is where the solver works.
* **`force` and `force!` return a force, `-∇U`, not the energy gradient**, and
  they do **not** divide by mass — that is the caller's job. Before this
  package was decoupled, the long-range entry point returned the gradient and
  the caller wrote `acceleration -= grad / mass`; it now returns `-grad` and the
  MD adapter writes `acceleration += force / mass`.
* **A coincident pair (`r == 0`) contributes exactly zero** rather than raising.
  It is reachable from an ordinary lattice initialisation (two sites whose x/y
  separation is an exact multiple of `Lx`/`Ly` at equal `z`), and the
  alternative is a `0/0` NaN in the force array that would poison every other
  particle.
* `force!` **fills** its buffer, it does not accumulate into it.

The accuracy references `direct_sum`, `soe_direct_sum` and `diff_direct_sum`
take the same `(plan, poses, charges)` arguments, so they work standalone too.
Note that they do not carry the `1 / (4π ϵ_0)` prefactor `SoEwald2D.energy`
applies, so they agree with it only at `ϵ_0 = 1`; that convention predates the
decoupling and is preserved.

## MD usage via ExTinyMD

`using ExTinyMD` loads the `SoEwald2DExTinyMDExt` extension, which supplies
`ExTinyMD.AbstractInteraction` wrappers around the plans plus the
`ExTinyMD.energy` and `ExTinyMD.update_acceleration!` methods for them.

```julia
using ExTinyMD, SoEwald2D

n_atoms = 100
L = 20.0
boundary = ExTinyMD.Q2dBoundary(L, L, L)

atoms = Vector{Atom{Float64}}()
for i in 1:n_atoms ÷ 2
    push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
end
for i in n_atoms ÷ 2 + 1 : n_atoms
    push!(atoms, Atom(type = 2, mass = 1.0, charge = -1.0))
end

info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.5, L - 0.5), boundary; min_r = 1.0, temp = 1.0)

ϵ_0 = 1.0
α = 0.3
s = 2.0
r_c = s / α            # 6.667 < min(Lx, Ly) / 2 = 10.0
k_c = 2 * s * α

no_finder = NoNeighborFinder()
celllist = CellList3D(info, r_c, boundary, 100)
interaction_short = SoEwald2DShortInteraction(ϵ_0, (L, L, L), s, α, n_atoms, r_c)
interaction_long  = SoEwald2DLongInteraction(ϵ_0, (L, L, L), s, α, n_atoms, k_c, SoePara())

interactions = [
    (LennardJones(), CellList3D(info, 4.5, boundary, 100)),
    (SubLennardJones(0.0, L; cutoff = 0.5, σ = 0.5), SubNeighborFinder(1.0, info, 0.0, L)),
    (interaction_short, celllist),
    (interaction_long, no_finder),
]
loggers = [TemperatureLogger(100, output = false), TrajectoryLogger(step = 100, output = false)]
simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

sys = MDSys(
    n_atoms = n_atoms,
    atoms = atoms,
    boundary = boundary,
    interactions = interactions,
    loggers = loggers,
    simulator = simulator,
)

simulate!(simulator, sys, info, 100)
```

The short-range wrapper needs a **3-D** neighbour finder — `CellList3D` or
`CellListDir3D` — because its kernel is `erfc(α·r)/r` in the full
three-dimensional separation. `NoNeighborFinder()` is also accepted and makes
the plan fall back to its own `O(n²)` pair loop. Anything else, including the
quasi-2D finders `CellListQ2D`/`CellListDirQ2D`, raises an `ArgumentError`
naming what to use instead: a mismatched finder would otherwise be silent. The
long-range wrapper ignores its finder entirely, so `NoNeighborFinder()` is the
natural choice there.

## Breaking changes in 0.2.0

Moving `ExTinyMD` to `[weakdeps]` is a breaking change, so 0.1.5 → 0.2.0.

* **`SoEwald2DShortInteraction` and `SoEwald2DLongInteraction` are now
  functions, not types.** They used to be `struct`s in `SoEwald2D`; a struct's
  supertype is fixed where the struct is defined and no extension can retrofit
  one, so the `ExTinyMD.AbstractInteraction` wrappers must live in the
  extension's own namespace. The two exported names are now *dispatcher
  functions* that forward to the real constructors there, which means:

  * **Construction works unchanged.** `SoEwald2DShortInteraction(ϵ_0, L, s, α,
    n_atoms, r_c)` returns exactly the wrapper it always did, and it still
    `isa ExTinyMD.AbstractInteraction`.
  * **Type-position uses do not work.** `x isa SoEwald2DShortInteraction`, an
    `::SoEwald2DLongInteraction` annotation, a
    `Vector{SoEwald2DShortInteraction}` element type, and dispatching a method
    on one all now raise a `TypeError`.

  If you need the type itself, reach into the extension module:

  ```julia
  using SoEwald2D, ExTinyMD
  ext = Base.get_extension(SoEwald2D, :SoEwald2DExTinyMDExt)
  x isa ext.SoEwald2DShortInteraction        # works
  ```

  Calling either name without `using ExTinyMD` raises an informative error
  rather than an `UndefVarError`.

* **`SoEwald2D_El`, `SoEwald2D_Es`, `SoEwald2D_Fl!` and `SoEwald2D_Fs!` are
  gone.** They were the `(interaction, sys, info)` adapter entry points. Use
  `SoEwald2D.energy`/`force`/`force!` on a plan for standalone work, or
  `ExTinyMD.energy`/`ExTinyMD.update_acceleration!` on a wrapper inside an
  `MDSys`.
* **`revise_interaction!` is gone** — gathering from `sys`/`info` is the
  extension's job now.
* **`SoEwald2D_init` is gone.** It referenced an `n_atoms` that was not one of
  its parameters, so every call raised an unconditional `UndefVarError`;
  nothing but the export list mentioned it.
* **`direct_sum`, `soe_direct_sum` and `diff_direct_sum`** take
  `(plan, poses, charges)` instead of `(interaction, sys, info)`.
* **The short-range `update_acceleration!` now divides by mass.** The old
  `SoEwald2D_Fs!` added a force straight into `acceleration` with no mass
  division, which every test in the old suite hid by using `mass = 1.0`. The
  long-range path always did divide. Results change only for non-unit masses,
  and only in the direction of being correct.
* **`Base.real(::ExTinyMD.Point)` is gone.** It was type piracy on a type this
  package does not own, and it was unused.

## Installing from source / registration status

`SoEwald2D` 0.2.0 needs `ExTinyMD` 0.3, which is in the General registry, and
its **test** environment additionally needs `QuasiEwald` 0.3 as an
independent-implementation reference. `QuasiEwald` 0.3 is not registered yet,
so `test/Project.toml` carries a `[sources]` override pinning it to its
`decouple-extinymd` branch. General's automerge rejects any `Project.toml`
containing `[sources]`, so removing that pin gates the next release; the item
is tracked for all the affected packages together in
`ExTinyMD.jl/docs/superpowers/specs/2026-09-15-downstream-decoupling-design.md`,
§6b. Running the package itself needs no override.
