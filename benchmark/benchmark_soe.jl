using BenchmarkTools, SoEwald2D, Plots, LaTeXStrings, SpecialFunctions

begin
    Z = 0.0:0.01:100.0
    α = 2.0
    k = 1.5
    exact = exp.(k .* Z) .* erfc.(α .* Z .+ k / (2 * α))
    soepara = SoePara()
    soexpara = SoExPara()
    soerfc_result = [exp(k * z) * soerfc(α * z + k / (2 * α), soepara) for z in Z]
    soexp_result = [exp_mul_erfc(z, α, k, soexpara) for z in Z]
    plot(Z, log10.(abs.(exact .- soerfc_result)), label = "SoErfc")
    plot!(Z, log10.(abs.(exact .- soexp_result)), label = "SoExp")
    title!(L"$\exp{(k z)} \mathrm{erfc}(\alpha z + \frac{k}{2 \alpha})$")
    xlabel!(L"z")
    ylabel!(L"$\log{\mathrm{(Error)}}$")
    savefig("exp_mul_erfc_error.pdf")
end