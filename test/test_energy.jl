@testset "compare energy by SoEwald2D and DirectEwald2D" begin
    n_atoms = 1000
    L = (10.0, 10.0, 10.0)
    q = 2 .* rand(n_atoms) .- 1.0
    q = q .- sum(q) / n_atoms
    x = L[1] .* rand(n_atoms)
    y = L[2] .* rand(n_atoms)
    z = L[3] .* rand(n_atoms)

    para = SoEwald2DPara(L, rand(), rand(), n_atoms)

    iterpara = IterPara(n_atoms)
    soepara = SoePara()

    U_soe = [0.0]
    U_dir = [0.0]

    energy_sum!(q, x, y, z, para, soepara, iterpara, U_soe)
    direct_sum!(q, x, y, z, para, U_dir)

    @show U_soe[1], U_dir[1]
    @test U_soe[1] ≈ U_dir[1]
end