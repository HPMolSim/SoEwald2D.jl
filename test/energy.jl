@testset "compare energy by SoEwald2D and DirectEwald2D" begin
    n_atoms = 200
    L = 10.0
    boundary = ExTinyMD.Q2dBoundary(L, L, 10.0)

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms ÷ 2
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    for i in n_atoms ÷ 2 + 1: n_atoms
        push!(atoms, Atom(type = 2, mass = 1.0, charge = - 1.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, 10.0), boundary; min_r = 1.0, temp = 1.0)

    para = SoEwald2DLongInteraction(1.0, (L, L, L), 1e-5, 1.0, n_atoms, sqrt(- 4 * log(1e-4)), SoePara())

    interactions = [(LennardJones(), CellListDir3D(info, 4.5, boundary, 100))]
    loggers = [TempartureLogger(100, output = false)]
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

    sys = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers,
        simulator = simulator
    )

    U_soe = SoEwald2D_El(para, sys, info)
    U_dir = direct_sum(para, sys, info)
    U_soe_dir = soe_direct_sum(para, sys, info, para.soepara)

    @test U_soe ≈ U_dir
    @test isapprox(U_dir, U_soe_dir, atol = 1e-12)
end

@testset "compare energy with ICM" begin
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

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)

    interactions = [(LennardJones(), CellListDir3D(info, 4.5, boundary, 100))]
    loggers = [TempartureLogger(100, output = false)]
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

    sys = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers,
        simulator = simulator
    )
        
    ϵ_0 = 1.0
    accuracy = 1e-6
    α = 0.5
    r_c = 9.9
    k_c = sqrt(-4 * α * log(accuracy))

    no_finder = NoNeighborFinder(n_atoms);
    celllist = CellList3D(info, r_c, boundary, 1);
    interaction_long = SoEwald2DLongInteraction(ϵ_0, (L, L, L), accuracy, α, n_atoms, k_c, SoePara());
    interaction_short = SoEwald2DShortInteraction(ϵ_0, (L, L, L), accuracy, α, n_atoms, r_c);
    
    Es = ExTinyMD.energy(interaction_short, celllist, sys, info)
    El = ExTinyMD.energy(interaction_long, no_finder, sys, info)

    N_real = 100
    N_img = 0
    ICM_sys = IcmSys((0.0, 0.0), (L, L, L), N_real, N_img)
    coords = [p_info.position for p_info in info.particle_info]
    charge = [atoms[p_info.id].charge for p_info in info.particle_info]
    ref_pos, ref_charge = IcmSysInit(ICM_sys, coords, charge)
    energy_icm = IcmEnergy(ICM_sys, coords, charge, ref_pos, ref_charge)

    @test isapprox(energy_icm, (Es + El), atol = 1e-2)
end