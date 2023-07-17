
@testset "compare ad_energy_sum with ad_direct_sum and diff_direct_sum" begin
    n_atoms = 1000
    q = 2 .* rand(n_atoms) .- 1.0
    q = q .- sum(q) / n_atoms

    x = 10.0 * rand(n_atoms)
    y = 10.0 * rand(n_atoms)
    z = 10.0 * rand(n_atoms)

    U = [0.0]
    F_x = zeros(n_atoms)
    F_y = zeros(n_atoms)
    F_z = zeros(n_atoms)

    para = SoEwald2DPara((10.0, 10.0, 10.0), 1.0, rand(), n_atoms);
    soepara = SoePara()
    iterpara = IterPara(n_atoms)
    adpara_dir = AdPara(Float64, n_atoms)
    adpara_soe = AdPara(Float64, n_atoms)

    direct_sum!(q, x, y, z, para, U)
    diff_direct_sum!(q, x, y, z, para, F_x, F_y, F_z)
    ad_direct_sum!(q, x, y, z, para, adpara_dir)
    ad_energy_sum!(q, x, y, z, para, soepara, iterpara, adpara_soe)

    @test abs(adpara_dir.U[1] - adpara_soe.U[1]) < 1e-6
    @test sum(abs.(adpara_dir.Fx .- adpara_soe.Fx))/n_atoms < 1e-6
    @test sum(abs.(adpara_dir.Fy .- adpara_soe.Fy))/n_atoms < 1e-6
    @test sum(abs.(adpara_dir.Fz .- adpara_soe.Fz))/n_atoms < 1e-6

    @test abs(adpara_dir.U[1] - U[1]) < 1e-6
    @test sum(abs.(adpara_dir.Fx ./ 2 .- F_x))/n_atoms < 1e-6
    @test sum(abs.(adpara_dir.Fy ./ 2 .- F_y))/n_atoms < 1e-6
    @test sum(abs.(adpara_dir.Fz ./ 2 .- F_z))/n_atoms < 1e-6
end