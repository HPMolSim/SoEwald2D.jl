# @testset "simulation with ExTinyMD.jl, rbm off" begin
#     n_atoms = 100
#     L = 100.0
#     boundary = ExTinyMD.Q2dBoundary(L, L, L)

#     atoms = Vector{Atom{Float64}}()
#     for i in 1:n_atoms÷2
#         push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
#     end

#     for i in n_atoms÷2 + 1 : n_atoms
#         push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
#     end

#     info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.5, L - 0.5), boundary; min_r = 1.0, temp = 1.0)

#     ϵ_0 = 1.0 / 3.5
#     α = 0.1
#     s = 1.0
#     r_c = s / α
#     k_c = 2 * s * α

#     no_finder = NoNeighborFinder(n_atoms);
#     celllist = CellList3D(info, r_c, boundary, 1);
#     interaction_long = SoEwald2DLongInteraction(ϵ_0, (L, L, L), s, α, n_atoms, k_c, SoePara());
#     interaction_short = SoEwald2DShortInteraction(ϵ_0, (L, L, L), s, α, n_atoms, r_c);

#     interactions = [
#         (LennardJones(), CellList3D(info, 4.5, boundary, 100)),
#         (SubLennardJones(0.0, L; cutoff = 0.5, σ = 0.5), SubNeighborFinder(1.0, info, 0.0, L)),
#         (interaction_short, celllist),
#         (interaction_long, no_finder)
#         ]
#     loggers = [TempartureLogger(100, output = false), TrajectionLogger(step = 100, output = false)]
#     simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

#     sys = MDSys(
#         n_atoms = n_atoms,
#         atoms = atoms,
#         boundary = boundary,
#         interactions = interactions,
#         loggers = loggers,
#         simulator = simulator
#     )
    
#     simulate!(simulator, sys, info, 10)

#     @test info.running_step == 10
# end

@testset "simulation with ExTinyMD.jl, rbm on" begin
    n_atoms = 100
    L = 100.0
    boundary = ExTinyMD.Q2dBoundary(L, L, L)

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms÷2
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    for i in n_atoms÷2 + 1 : n_atoms
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.5, L - 0.5), boundary; min_r = 1.0, temp = 1.0)

    ϵ_0 = 1.0 / 3.5
    α = 0.1
    s = 3.0
    r_c = s / α
    k_c = 2 * s * α

    no_finder = NoNeighborFinder(n_atoms);
    celllist = CellList3D(info, r_c, boundary, 1);
    interaction_long = SoEwald2DLongInteraction(ϵ_0, (L, L, L), s, α, n_atoms, k_c, SoePara(), rbm = true, rbm_p = 30);
    interaction_short = SoEwald2DShortInteraction(ϵ_0, (L, L, L), s, α, n_atoms, r_c);

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
    
    simulate!(simulator, sys, info, 10)

    @test info.running_step == 10 
end