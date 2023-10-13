# SoEwald2D

[![Build Status](https://github.com/ArrogantGao/SoEwald2D.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ArrogantGao/SoEwald2D.jl/actions/workflows/CI.yml?query=branch%3Amain)

`SoEwald2D.jl` is a package based on `Julia` for rapid calculation of electrostatic interaction for particles quasi-2D charged systems.
It combined the SOE method and the Ewald2D method so that the complexity of its pairwise summation in the k space is of $O(N)$, much better than the original one, which is of $O(N^2)$.

This package can work as a force field of the MD package `ExTinyMD.jl`.

## Getting Started

Install the package in Julia REPL simply by
```julia
pkg> add ExTinyMD, SoEwald2D
```
And here is an example about using it as the force field of `ExTinyMD.jl`:
```julia
using ExTinyMD, SoEwald2D

begin
    n_atoms = 100
    L = 20.0
    boundary = ExTinyMD.Q2dBoundary(L, L, L)
    
    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms÷2
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end
    
    for i in n_atoms÷2 + 1 : n_atoms
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
    end
    
    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.5, L - 0.5), boundary; min_r = 1.0, temp = 1.0)
    
    ϵ_0 = 1.0
    accuracy = 1e-4
    α = 0.5
    r_c = 5.0
    k_c = sqrt(- 4 * α * log(accuracy))
    
    no_finder = NoNeighborFinder(n_atoms);
    celllist = CellList3D(info, r_c, boundary, 100);
    interaction_long = SoEwald2DLongInteraction(ϵ_0, (L, L, L), accuracy, α, n_atoms, k_c, SoePara());
    interaction_short = SoEwald2DShortInteraction(ϵ_0, (L, L, L), accuracy, α, n_atoms, r_c);
    
    interactions = [
        (LennardJones(), CellList3D(info, 4.5, boundary, 100)),
        (SubLennardJones(0.0, L; cutoff = 0.5, σ = 0.5), SubNeighborFinder(1.0, info, 0.0, L)),
        (interaction_short, celllist),
        (interaction_long, no_finder)
        ]
    loggers = [TempartureLogger(100, output = false), TrajectionLogger(step = 100, output = false)]
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))
    
    sys = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers,
        simulator = simulator
    )
    
    simulate!(simulator, sys, info, 100)
end
```
