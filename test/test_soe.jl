@testset begin "SOE and exp_mul_erfc by SOE"
    soepara = SoePara()
    k_array = 0.01:0.01:100.0
    α = 10.0 * rand()
    z = 10.0 * rand()

    count_2 = 0
    count_3 = 0
    for k in k_array
        result_1 = zero(ComplexF64)
        result_2 = zero(ComplexF64)
        result_3 = zero(ComplexF64)

        for (s, w) in soepara.sw
            result_1 += s * w * exp(-k^2/(4 * α^2)) * α * exp( - s * α * z) / (s * α + k)
            if k/(2α) - α * z > 0
                result_2 += w * exp( - k * z) * exp( - s * k / (2α)) * exp(s * α * z)
            else
                result_3 += w * exp( - k * z) * exp(s * k / (2α)) * exp(- s * α * z)
            end
        end

        @test abs(result_1 - exp(k * z) * erfc(k/(2α) + α * z)) < 1e-12
        if k/(2α) - α * z > 0
            @test abs(result_2 - exp(-k * z) * erfc(k/(2α) - α * z)) < 1e-12
        elseif k/(2α) - α * z < 0
            @test abs(result_3 - exp(-k * z) * erfc(α * z - k/(2α))) < 1e-12
        end
    end
end