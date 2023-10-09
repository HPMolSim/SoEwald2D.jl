@testset "m list test" begin
    n_atoms = 200
    iterpara = IterPara(n_atoms)
    z = 10 .* rand(n_atoms)
    d = rand() + 1.0
    update_iterpara_z!(iterpara, z)
    update_iterpara_m!(iterpara, z, d)
    z_sort = sort(z)
    for i in 1:n_atoms
        m_i = iterpara.m_list[i]
        if m_i == 1
            @test z[iterpara.z_list[i]] - z[iterpara.z_list[m_i]] < d
        else
            @test z[iterpara.z_list[i]] - z[iterpara.z_list[m_i]] < d && z[iterpara.z_list[i]] - z[iterpara.z_list[m_i - 1]] > d
        end
    end
end

@testset "difference between direct_iter and iter" begin
    n_atoms = 1000
    L = (10.0, 10.0, 10.0)
    q = 2 .* rand(n_atoms) .- 1.0
    q = q .- sum(q) / n_atoms
    x = L[1] .* rand(n_atoms)
    y = L[2] .* rand(n_atoms)
    z = L[3] .* rand(n_atoms)

    para = SoEwald2DPara(L, 3.0, rand(), n_atoms)

    iterpara = IterPara(n_atoms)
    soepara = SoePara()

    for k in 0.01:0.5:10.0
        K = (k, 0.0, k)
        α = para.α
        update_iterpara_z!(iterpara, z)
        update_iterpara_m!(iterpara, z, k / (2 * α^2))

        sum_k = zero(ComplexF64)

        for (s, w) in soepara.sw
            D_dir = direct_iterpara_D(q, x, y, z, para, s, K)
            update_iterpara_D!(iterpara, q, x, y, z, para, s, K)
            @test D_dir ≈ iterpara.D
        end
    end
end