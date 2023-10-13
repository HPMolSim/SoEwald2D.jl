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