@testset begin "SOE and exp_mul_erfc by SOE"
    soepara4 = SoePara4()
    soepara8 = SoePara8()
    soepara16 = SoePara()
    α = 1.0
    z = 1.0
    
    for k in 0.01:0.1:10.0
        for (soepara, atol) in zip([soepara4, soepara8, soepara16], [1e-3, 1e-7, 1e-13])
            @test isapprox(soexp(k, soepara), exp(-k^2), atol=atol)
            @test isapprox(soexp(-k, soepara), exp(-k^2), atol=atol)
            @test isapprox(soerf(k, soepara), erf(k), atol=atol)
            @test isapprox(soerf(-k, soepara), erf(-k), atol=atol)
            @test isapprox(soerfc(k, soepara), erfc(k), atol=atol)
            @test isapprox(soerfc(-k, soepara), erfc(-k), atol=atol)
            @test isapprox(soexp_mul_erfc(z, α, k, soepara), exp(k * z) * erfc(k / (2 * α) + α * z), atol=atol)
        end
    end
end