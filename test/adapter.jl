@testset "ExTinyMD adapter (Task 4): the wrapper under load" begin
    # The two tests here catch different things and neither substitutes for the
    # other. The permutation test catches id/slot and mass-division faults that
    # a trajectory cannot see; the conservation test catches force faults that
    # a per-call comparison cannot, because it is the only place the integrator
    # actually uses the forces.

    function _charged_system(n, L, Lz; temp = 1.0)
        boundary = Q2dBoundary(L, L, Lz)
        atoms = Atom{Float64}[]
        for i in 1:(n ÷ 2)
            push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
        end
        for i in (n ÷ 2 + 1):n
            push!(atoms, Atom(type = 2, mass = 1.0, charge = -1.0))
        end
        info = SimulationInfo(n, atoms, (0.0, L, 0.0, L, 1.0, Lz - 1.0), boundary; min_r = 1.0, temp = temp)
        info.running_step = 1
        return boundary, atoms, info
    end

    # Every id gets a distinct mass and a distinct charge magnitude, so that an
    # id/slot mix-up or a missing mass division cannot cancel out. Velocities
    # were drawn against the uniform-mass atoms above; that changes the initial
    # condition only, not the conservation law being tested.
    _nonuniform(atoms, n) = [Atom(type = a.type, mass = 0.5 + 0.15 * i,
                                  charge = (i <= n ÷ 2 ? 1.0 : -1.0) * (1 + 0.02 * i))
                             for (i, a) in enumerate(atoms)]

    _kinetic(atoms, info) = sum(0.5 * atoms[p.id].mass *
                                (p.velocity[1]^2 + p.velocity[2]^2 + p.velocity[3]^2)
                                for p in info.particle_info)

    @testset "adapter is correct when slot order differs from id order" begin
        # In stock ExTinyMD `particle_info[i].id == i`, so a gather that
        # confuses slot with id looks correct forever. Permute the mapping so
        # the two differ, following ExTinyMD's own test_adapter.jl pattern:
        # reverse the slot order and keep ids attached to their particles via
        # info.id_dict.
        Random.seed!(20260939)
        n, L, Lz = 12, 10.0, 8.0
        boundary, atoms, info = _charged_system(n, L, Lz)
        atoms = _nonuniform(atoms, n)

        reverse!(info.particle_info)
        for i in eachindex(info.particle_info)
            info.id_dict[info.particle_info[i].id] = i
        end
        @test info.particle_info[1].id != 1   # the mapping really is permuted
        @test length(unique(a.mass for a in atoms)) == n     # and every mass
        @test length(unique(a.charge for a in atoms)) == n   # and charge differs

        ϵ_0 = 1.0
        α = 0.4
        s = 2.0
        # s / α = 5.0 is NOT < min(Lx, Ly) / 2 = 5.0 at this box size, and
        # SoEwald2DShortPlan rejects that, so the cutoff is set independently.
        r_c = 4.0           # r_c = 4.0 < min(Lx, Ly) / 2 = 5.0
        k_c = 2 * s * α

        intershort = SoEwald2DShortInteraction(ϵ_0, (L, L, Lz), s, α, n, r_c)
        short_finder = CellList3D(info, r_c, boundary, 1)
        interlong = SoEwald2DLongInteraction(ϵ_0, (L, L, Lz), s, α, n, k_c, SoePara(); parallel = false)
        long_finder = NoNeighborFinder()

        sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                    interactions = [(intershort, short_finder), (interlong, long_finder)],
                    loggers = [TemperatureLogger(10^9; output = false)],
                    simulator = VerletProcess(dt = 0.001))

        # The reference: gathered by hand in SLOT order, with charge looked up
        # by ID. If the adapter did it the other way round these would differ.
        poses = [SVector(p.position[1], p.position[2], p.position[3]) for p in info.particle_info]
        charges = [atoms[p.id].charge for p in info.particle_info]

        Es_ref = SoEwald2D.energy(intershort.plan, poses, charges; neighbor_list = short_finder.neighbor_list)
        El_ref = SoEwald2D.energy(interlong.plan, poses, charges)

        # rtol 1e-12: every fault this is here to catch (id/slot mix-up,
        # missing mass division, factor of two) is an O(1) relative error,
        # twelve orders of magnitude above the bound. Measured deviation with
        # the adapter correct: 0.
        @test isapprox(ExTinyMD.energy(intershort, short_finder, sys, info), Es_ref, rtol = 1e-12)
        @test isapprox(ExTinyMD.energy(interlong, long_finder, sys, info), El_ref, rtol = 1e-12)

        Fs_ref = SoEwald2D.force(intershort.plan, poses, charges; neighbor_list = short_finder.neighbor_list)
        Fl_ref = SoEwald2D.force(interlong.plan, poses, charges)

        for p in info.particle_info
            p.acceleration = Point(0.0, 0.0, 0.0)
        end
        ExTinyMD.update_acceleration!(intershort, short_finder, sys, info)
        ExTinyMD.update_acceleration!(interlong, long_finder, sys, info)

        for (slot, p) in enumerate(info.particle_info)
            m = atoms[p.id].mass
            for d in 1:3
                @test isapprox(p.acceleration[d],
                               Fs_ref[slot][d] / m + Fl_ref[slot][d] / m, rtol = 1e-12)
            end
        end

        # update_acceleration! must ACCUMULATE, not overwrite: running it again
        # on top of what is already there has to double the result.
        doubled = [p.acceleration for p in info.particle_info]
        ExTinyMD.update_acceleration!(intershort, short_finder, sys, info)
        ExTinyMD.update_acceleration!(interlong, long_finder, sys, info)
        for (slot, p) in enumerate(info.particle_info), d in 1:3
            @test isapprox(p.acceleration[d], 2 * doubled[slot][d], rtol = 1e-12)
        end
    end

    @testset "adapter conserves KE + E_elec inside simulate!" begin
        # WHY THIS ASSERTS ON KE + E_elec AND NOT ON E_elec ALONE.
        #
        # A bound on the electrostatic energy drift cannot catch a force error.
        # Measured over the same eight seeds used to set the bound below, with
        # an exact 2x electrostatic force fault (both (interaction, finder)
        # pairs listed twice in sys.interactions):
        #
        #   |ΔE_elec|  correct 2.17e-2 ... 1.41e-1
        #   |ΔE_elec|  2x fault 1.31e-2 ... 1.68e-1
        #
        # The two ranges OVERLAP -- on one seed the faulted run drifted *less*
        # than the correct one. E_elec after N steps is set by thermal motion,
        # not by whether the force is right, so no bound on it can separate the
        # two cases.
        #
        # The quantity that does respond is the conserved one. `VerletProcess`
        # with the default `NoThermoStat` is symplectic and electrostatics is
        # the only interaction in `sys`, so KE + E_elec is conserved to O(dt^2).
        # Double the force and the integrator conserves KE + 2*E_elec instead,
        # so KE + E_elec drifts by the electrostatic work -- first order, not
        # second.
        #
        #   |Δ(KE + E_elec)|  correct   7.59e-5 ... 1.18e-3   (8 seeds)
        #   |Δ(KE + E_elec)|  2x fault  1.54e-2 ... 1.68e-1   (same 8 seeds)
        #
        # Bound 5e-3: 4.2x above the worst correct drift over the eight seeds
        # (and 66x above this seed's own 7.59e-5), 3.1x below the smallest
        # faulted drift over the same eight (and 31x below this seed's own
        # 1.57e-1). Verified by sabotage to fail on the 2x fault and pass
        # without it.
        #
        # temp = 0.05 rather than 1.0 keeps the slab from flying apart over 500
        # steps at this dt, there being no thermostat to hold it.
        Random.seed!(20260937)
        n, L, Lz = 20, 12.0, 10.0
        boundary, atoms, info = _charged_system(n, L, Lz; temp = 0.05)
        atoms = _nonuniform(atoms, n)

        ϵ_0 = 1.0
        α = 0.4
        s = 2.0
        r_c = s / α         # r_c = 5.0 < min(Lx, Ly) / 2 = 6.0
        k_c = 2 * s * α
        @test r_c < min(L, L) / 2

        intershort = SoEwald2DShortInteraction(ϵ_0, (L, L, Lz), s, α, n, r_c)
        short_finder = CellList3D(info, r_c, boundary, 1)
        interlong = SoEwald2DLongInteraction(ϵ_0, (L, L, Lz), s, α, n, k_c, SoePara(); parallel = false)
        long_finder = NoNeighborFinder()

        sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                    interactions = [(intershort, short_finder), (interlong, long_finder)],
                    loggers = [TemperatureLogger(10^9; output = false)],
                    simulator = VerletProcess(dt = 5e-3))

        Eel0 = ExTinyMD.energy(intershort, short_finder, sys, info) +
               ExTinyMD.energy(interlong, long_finder, sys, info)
        K0 = _kinetic(atoms, info)
        simulate!(sys.simulator, sys, info, 500)
        Eel1 = ExTinyMD.energy(intershort, short_finder, sys, info) +
               ExTinyMD.energy(interlong, long_finder, sys, info)
        K1 = _kinetic(atoms, info)

        @test isfinite(Eel1)
        @test isfinite(K1)
        # The load-bearing assertion: the forces the integrator used must be
        # the gradient of the energy these same wrappers report. A sign error,
        # a scaling error, a dropped or double-counted mass division, and a
        # dropped or double-counted pair all break that equality and show up
        # here.
        @test abs((K1 + Eel1) - (K0 + Eel0)) < 5e-3
        # The electrostatic energy itself must at least stay in range. Smoke
        # check only, for the reason given above: measured up to 1.41e-1 over
        # the eight seeds, so 1.0 is a 7x margin and it cannot distinguish a
        # correct run from a faulted one.
        @test abs(Eel1 - Eel0) < 1.0
    end

    @testset "the short-range wrapper refuses non-3-D neighbour finders" begin
        # An untyped `_finder_list(f) = f.neighbor_list` accepts CellListQ2D /
        # CellListDirQ2D happily, because they have a field of that name -- and
        # then feeds a sum whose cutoff is on the 3-D separation a candidate
        # set chosen by in-plane distance. Before the decoupling those were a
        # MethodError; the fallback method keeps the failure loud.
        Random.seed!(20260940)
        n, L, Lz = 12, 10.0, 8.0
        boundary, atoms, info = _charged_system(n, L, Lz)

        ϵ_0, α, s = 1.0, 0.4, 2.0
        r_c = 4.0           # r_c = 4.0 < min(Lx, Ly) / 2 = 5.0

        intershort = SoEwald2DShortInteraction(ϵ_0, (L, L, Lz), s, α, n, r_c)
        sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                    interactions = [(intershort, CellList3D(info, r_c, boundary, 1))],
                    loggers = [TemperatureLogger(10^9; output = false)],
                    simulator = VerletProcess(dt = 0.001))

        for bad in (CellListQ2D(info, r_c, boundary, 1), CellListDirQ2D(info, r_c, boundary, 1))
            @test_throws ArgumentError ExTinyMD.energy(intershort, bad, sys, info)
            @test_throws ArgumentError ExTinyMD.update_acceleration!(intershort, bad, sys, info)
        end

        # The supported finders all still work, and -- since the plan always
        # recomputes the true separation from the candidate list -- the two 3-D
        # finders and the O(n^2) fallback must agree. rtol 1e-12: the only
        # difference between them is summation order.
        E_cl = ExTinyMD.energy(intershort, CellList3D(info, r_c, boundary, 1), sys, info)
        E_dir = ExTinyMD.energy(intershort, CellListDir3D(info, r_c, boundary, 1), sys, info)
        E_none = ExTinyMD.energy(intershort, NoNeighborFinder(), sys, info)
        @test isapprox(E_cl, E_none, rtol = 1e-12)
        @test isapprox(E_dir, E_none, rtol = 1e-12)

        # And the message has to say what to use, not just that it failed.
        msg = try
            ExTinyMD.energy(intershort, CellListQ2D(info, r_c, boundary, 1), sys, info)
            ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("CellList3D", msg)
        @test occursin("CellListQ2D", msg)
    end

    @testset "the two preserved names are functions, not types" begin
        # Documented breaking change: `SoEwald2DShortInteraction` and
        # `SoEwald2DLongInteraction` used to be structs and are now dispatcher
        # functions, so construction is unchanged but every type-position use
        # of the bare name raises a TypeError. This test is the executable form
        # of that note in the README -- if the pattern is ever changed back,
        # this is what says the documentation has gone stale.
        Random.seed!(20260941)
        n, L, Lz = 6, 10.0, 8.0
        boundary, atoms, info = _charged_system(n, L, Lz)
        ϵ_0, α, s = 1.0, 0.4, 2.0
        r_c = 4.0           # r_c = 4.0 < min(Lx, Ly) / 2 = 5.0

        inter = SoEwald2DShortInteraction(ϵ_0, (L, L, Lz), s, α, n, r_c)
        interl = SoEwald2DLongInteraction(ϵ_0, (L, L, Lz), s, α, n, 2 * s * α, SoePara())

        # Construction works exactly as before, and the result still subtypes
        # ExTinyMD's abstract type.
        @test inter isa ExTinyMD.AbstractInteraction
        @test interl isa ExTinyMD.AbstractInteraction
        # ... and it holds the framework-free plan.
        @test inter.plan isa SoEwald2DShortPlan
        @test interl.plan isa SoEwald2DLongPlan

        # But the exported names are not types.
        @test !(SoEwald2DShortInteraction isa Type)
        @test !(SoEwald2DLongInteraction isa Type)
        @test_throws TypeError (inter isa SoEwald2DShortInteraction)
        @test_throws TypeError (interl isa SoEwald2DLongInteraction)

        # The documented workaround.
        ext = Base.get_extension(SoEwald2D, :SoEwald2DExTinyMDExt)
        @test ext !== nothing
        @test inter isa ext.SoEwald2DShortInteraction
        @test interl isa ext.SoEwald2DLongInteraction

        # Constructing from a plan directly also works.
        @test SoEwald2DShortInteraction(SoEwald2DShortPlan(ϵ_0, (L, L, Lz), s, α, n, r_c)) isa
              ext.SoEwald2DShortInteraction
    end
end
