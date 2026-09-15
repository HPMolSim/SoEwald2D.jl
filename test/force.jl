@testset "compare the force_sum with diff_direct_sum" begin
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

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)

    interactions = [(LennardJones(), CellList3D(info, 4.5, boundary, 100))]
    loggers = [TemperatureLogger(100, output = false)]
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
    s = 1.0
    α = 0.1
    r_c = s / α       # r_c = 10.0 < min(Lx, Ly) / 2 = 50.0
    k_c = 2 * s * α

    plan = SoEwald2DLongPlan(ϵ_0, (L, L, L), s, α, n_atoms, k_c, SoePara());

    poses = [SVector(p_info.position[1], p_info.position[2], p_info.position[3]) for p_info in info.particle_info]
    charges = [atoms[p_info.id].charge for p_info in info.particle_info]

    # `diff_direct_sum` returns -grad U, the same sign as `SoEwald2D.force`.
    # Before the decoupling this comparison went through
    # `info.particle_info[i].acceleration` after `SoEwald2D_Fl!`, which was
    # -grad U / mass at the unit masses used here -- so the quantity being
    # compared is unchanged, only the route to it is.
    F = SoEwald2D.force(plan, poses, charges)
    sum_direct = diff_direct_sum(plan, poses, charges)

    for i in 1:n_atoms
        @test sum_direct[i][1] ≈ F[i][1]
        @test sum_direct[i][2] ≈ F[i][2]
        @test sum_direct[i][3] ≈ F[i][3]
    end
end

@testset "compare sort sum with ICM" begin
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

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)

    interactions = [(LennardJones(), CellListDir3D(info, 4.5, boundary, 100))]
    loggers = [TemperatureLogger(100, output = false)]
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

    sys = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers,
        simulator = simulator
    )
        
    ϵ_0 = 1.0 / 3.5
    α = 0.1
    s = 3.0
    r_c = s / α       # r_c = 30.0 < min(Lx, Ly) / 2 = 50.0
    k_c = 2 * s * α

    no_finder = NoNeighborFinder();
    celllist = CellList3D(info, r_c, boundary, 1);
    interaction_long = SoEwald2DLongInteraction(ϵ_0, (L, L, L), s, α, n_atoms, k_c, SoePara());
    interaction_short = SoEwald2DShortInteraction(ϵ_0, (L, L, L), s, α, n_atoms, r_c);

    ExTinyMD.update_acceleration!(interaction_short, celllist, sys, info)
    ExTinyMD.update_acceleration!(interaction_long, no_finder, sys, info)


    N_real = 200
    N_img = 0
    ICM_sys = IcmSys((0.0, 0.0), (L, L, L), N_real, N_img)
    coords = [p_info.position for p_info in info.particle_info]
    charge = [atoms[p_info.id].charge for p_info in info.particle_info]
    ref_pos, ref_charge = IcmSysInit(ICM_sys, coords, charge)
    force_icm = IcmForce(ICM_sys, coords, charge, ref_pos, ref_charge) ./ ϵ_0

    for i in 1:n_atoms
        # QuasiEwald 0.3 (decoupled) returns SVector{3} from IcmForce, while
        # `acceleration` is an ExTinyMD `Point`, so `dist2` no longer has a
        # method for the pair. Index both instead -- the one interface the two
        # types share.
        acc = info.particle_info[i].acceleration
        f = force_icm[i]
        error_i = sqrt((f[1] - acc[1])^2 + (f[2] - acc[2])^2 + (f[3] - acc[3])^2)
        @test error_i < 1e-3
    end
end

@testset "non-netural system" begin
    n_atoms = 100
    L = 100.0
    boundary = ExTinyMD.Q2dBoundary(L, L, L)

    atoms = Vector{Atom{Float64}}()
    for i in 1:n_atoms
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 1.0, temp = 1.0)

    interactions = [(LennardJones(), CellListDir3D(info, 4.5, boundary, 100))]
    loggers = [TemperatureLogger(100, output = false)]
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

    sys = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers,
        simulator = simulator
    )
        
    ϵ_0 = 1.0 / 3.5
    α = 0.1
    s = 3.0
    r_c = s / α       # r_c = 30.0 < min(Lx, Ly) / 2 = 50.0
    k_c = 2 * s * α

    no_finder = NoNeighborFinder();
    celllist = CellList3D(info, r_c, boundary, 1);
    interaction_long = SoEwald2DLongInteraction(ϵ_0, (L, L, L), s, α, n_atoms, k_c, SoePara());
    interaction_short = SoEwald2DShortInteraction(ϵ_0, (L, L, L), s, α, n_atoms, r_c);

    ExTinyMD.update_acceleration!(interaction_short, celllist, sys, info)
    ExTinyMD.update_acceleration!(interaction_long, no_finder, sys, info)


    N_real = 200
    N_img = 0
    ICM_sys = IcmSys((0.0, 0.0), (L, L, L), N_real, N_img)
    coords = [p_info.position for p_info in info.particle_info]
    charge = [atoms[p_info.id].charge for p_info in info.particle_info]
    ref_pos, ref_charge = IcmSysInit(ICM_sys, coords, charge)
    force_icm = IcmForce(ICM_sys, coords, charge, ref_pos, ref_charge) ./ ϵ_0

    for i in 1:n_atoms
        # QuasiEwald 0.3 (decoupled) returns SVector{3} from IcmForce, while
        # `acceleration` is an ExTinyMD `Point`, so `dist2` no longer has a
        # method for the pair. Index both instead -- the one interface the two
        # types share.
        acc = info.particle_info[i].acceleration
        f = force_icm[i]
        error_i = sqrt((f[1] - acc[1])^2 + (f[2] - acc[2])^2 + (f[3] - acc[3])^2)
        @test error_i < 1e-3
    end
end