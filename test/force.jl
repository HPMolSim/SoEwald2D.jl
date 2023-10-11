@testset "compare the force_sum with diff_direct_sum" begin
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

    interactions = [(LennardJones(), CellList3D(info, 4.5, boundary, 100))]
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
    accuracy = 1e-5
    α = 2.0
    r_c = 4.5
    k_c = sqrt(-4 * α * log(accuracy))

    interaction = SoEwald2DLongInteraction(ϵ_0, (L, L, L), accuracy, α, n_atoms, k_c, SoePara());

    revise_interaction!(interaction, sys, info)

    iterpara = interaction.iterpara
    soepara = interaction.soepara
    q = interaction.q
    x = interaction.x
    y = interaction.y
    z = interaction.z
    mass = interaction.mass

    update_iterpara_z!(iterpara, z)

    # @testset "k_0" begin
    #     sum_sort = [Point(zero(ComplexF64), zero(ComplexF64), zero(ComplexF64)) for _=1:interaction.n_atoms]
    #     sum_direct = [Point(zero(Float64), zero(Float64), zero(Float64)) for _=1:interaction.n_atoms]

    #     force_sum_k0!(q, z, interaction, soepara, iterpara, sum_sort)
    #     diff_direct_sum_k0!(q, z, interaction, sum_direct)

    #     for i in 1:n_atoms
    #         @test real(sum_sort[i][3]) ≈ sum_direct[i][3]
    #     end
    # end

    @testset "k, xy" begin
        K = (0.3, 0.4, 0.5)

        sum_sort = [Point(zero(ComplexF64), zero(ComplexF64), zero(ComplexF64)) for _=1:interaction.n_atoms]
        sum_direct = [Point(zero(Float64), zero(Float64), zero(Float64)) for _=1:interaction.n_atoms]

        force_sum_k!(K, q, x, y, z, interaction, soepara, iterpara, sum_sort)
        diff_direct_sum_k!(K, q, x, y, z, interaction, sum_direct)

        for i in 1:n_atoms
            if real(sum_sort[i][1]) ≈ sum_direct[i][1]
                @show i, z[i], real(sum_sort[i]), sum_direct[i]
            end
            # @test real(sum_sort[i][1]) ≈ sum_direct[i][1]
            # @test real(sum_sort[i][2]) ≈ sum_direct[i][2]
        end
    end
end