@testset begin "SOE and exp_mul_erfc by SOE"
    soepara4 = SoePara4()
    soepara8 = SoePara8()
    soepara16 = SoePara()
    
    for k in 0.01:0.1:10.0
        @test isapprox(soerfc(k, soepara4), erfc(k), atol=1e-3)
        @test isapprox(soerfc(k, soepara8), erfc(k), atol=1e-7)
        @test isapprox(soerfc(k, soepara16), erfc(k), atol=1e-14)
    end
end