using BenchmarkTools, SoEwald2D, Plots, LaTeXStrings, SpecialFunctions

# benchmark time cost
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

#benchmark error of the summation
begin
    n_atoms = 100
    q = 2 .* rand(n_atoms) .- 1.0
    q = q .- sum(q) / n_atoms

    x = 10.0 * rand(n_atoms)
    y = 10.0 * rand(n_atoms)
    z = 10.0 * rand(n_atoms)

    U_soe = [0.0]
    U = [0.0]
    F_x = zeros(n_atoms)
    F_y = zeros(n_atoms)
    F_z = zeros(n_atoms)

    para = SoEwald2DPara((10.0, 10.0, 10.0), 1.0, rand(), n_atoms);
    soepara = SoePara()
    iterpara = IterPara(n_atoms)
    adpara_dir = AdPara(Float64, n_atoms)
    adpara_soe = AdPara(Float64, n_atoms)

    energy_sum!(q, x, y, z, para, soepara, iterpara, U_soe)
    direct_sum!(q, x, y, z, para, U)
    diff_direct_sum!(q, x, y, z, para, F_x, F_y, F_z)
    # ad_direct_sum!(q, x, y, z, para, adpara_dir)
    ad_energy_sum!(q, x, y, z, para, soepara, iterpara, adpara_soe)

    error_energy = abs(U_soe[1] - U[1])
    error_Force = sqrt(maximum((abs.(adpara_soe.Fx ./ 2 .- F_x)).^2 .+ (abs.(adpara_soe.Fy ./ 2 .- F_y)).^2 .+ (abs.(adpara_soe.Fz ./ 2 .- F_z)).^2))
    @show error_energy, error_Force
end

begin
    n_atoms = 100
    q = 2 .* rand(n_atoms) .- 1.0
    q = q .- sum(q) / n_atoms

    x = 10.0 * rand(n_atoms)
    y = 10.0 * rand(n_atoms)
    z = 10.0 * rand(n_atoms)

    para = SoEwald2DPara((10.0, 10.0, 10.0), 1.0, 1.0, n_atoms)
    soepara = SoePara()
    iterpara = IterPara(n_atoms)

    update_iterpara_z!(iterpara, z)

    # this will be used to benchmark the accuracy of the method
    k_array = 0.1:0.01:10.0
    Error_k = Vector{Float64}()
    for k in k_array
        K = (k, 0.0, k)
        sum_dir = direct_sum_k(K, q, x, y, z, para)
        sum_soe = energy_sum_k(K, q, x, y, z, para, soepara, iterpara)
        error = abs(sum_dir - sum_soe)
        push!(Error_k, error)
    end
    plot(k_array, log10.(Error_k))
end

begin
    n_atoms = 20
    L = (100.0, 100.0, 100.0)
    q = 2 .* rand(n_atoms) .- 1.0
    q = q .- sum(q) / n_atoms
    x = L[1] .* rand(n_atoms)
    y = L[2] .* rand(n_atoms)
    z = L[3] .* rand(n_atoms)

    para = SoEwald2DPara(L, rand(), rand(), n_atoms)

    iterpara = IterPara(n_atoms)
    soepara = SoePara()

    K = (3.0, 4.0, 5.0)
    soe_A = energy_sum_k_S1(K, q, x, y, z, para, soepara, iterpara)
    dir_soe_A, dir_A = sum_A(q, x, y, z, para, soepara, K)
    @show soe_A, dir_soe_A, dir_A
    @show soe_A - dir_A, dir_soe_A - dir_A
end

begin
    n_atoms = 100
    q = 2 .* rand(n_atoms) .- 1.0
    q = q .- sum(q) / n_atoms

    x = 10.0 * rand(n_atoms)
    y = 10.0 * rand(n_atoms)
    z = 10.0 * rand(n_atoms)

    para = SoEwald2DPara((10.0, 10.0, 10.0), 1.0, 1.0, n_atoms)
    soepara = SoePara()
    iterpara = IterPara(n_atoms)

    update_iterpara_z!(iterpara, z)

    # this will be used to benchmark the accuracy of the method
    k_array = 0.1:0.01:10.0
    Error_k_A = Vector{Float64}()
    Error_k_B = Vector{Float64}()
    Error_k_C = Vector{Float64}()
    Error_k_D = Vector{Float64}()
    Error_k = Vector{Float64}()
    for k in k_array
        K = (k, 0.0, k)
        ea, eb, ec, ed = abs.(direct_sum_k_ABCD(K, q, x, y, z, para) .- energy_sum_k_ABCD(K, q, x, y, z, para, soepara, iterpara))
        push!(Error_k_A, ea)
        push!(Error_k_B, eb)
        push!(Error_k_C, ec)
        push!(Error_k_D, ed)

        sum_dir = direct_sum_k(K, q, x, y, z, para)
        sum_soe = energy_sum_k(K, q, x, y, z, para, soepara, iterpara)
        error = abs(sum_dir - sum_soe)
        push!(Error_k, error)
    end
    plot(10.0 .* k_array, log10.(Error_k_A), label = "error of sum A")
    plot!(10.0 .* k_array, log10.(Error_k_B), label = "error of sum B")
    plot!(10.0 .* k_array, log10.(Error_k_C), label = "error of sum C")
    plot!(10.0 .* k_array, log10.(Error_k_D), label = "error of sum D")
    plot!(xlabel=L"$k L_z$", ylabel = L"$\log{(error)}$")
    # plot!(k_array, log10.(Error_k))
end

begin
    n_atoms = 100
    q = 2 .* rand(n_atoms) .- 1.0
    q = q .- sum(q) / n_atoms

    x = 10.0 * rand(n_atoms)
    y = 10.0 * rand(n_atoms)
    z = 10.0 * rand(n_atoms)

    para = SoEwald2DPara((10.0, 10.0, 10.0), 1.0, 1.0, n_atoms)
    soepara = SoePara()
    iterpara = IterPara(n_atoms)

    update_iterpara_z!(iterpara, z)

    # this will be used to benchmark the accuracy of the method
    k_array = 0.1:0.01:10.0
    Error_k_A = Vector{Float64}()
    Error_k_B = Vector{Float64}()
    Error_k_C = Vector{Float64}()
    Error_k_D = Vector{Float64}()
    Error_k = Vector{Float64}()
    for k in k_array
        K = (k, 0.0, k)
        da, db, dc, dd = direct_sum_k_ABCD(K, q, x, y, z, para)
        sa, sb, sc, sd = direct_sum_k_ABCD_soerfc(K, q, x, y, z, para, soepara)
        did, did2 = direct_iter_sum_D(K, q, x, y, z, para, soepara, iterpara)

        push!(Error_k_A, abs(dd - sd))
        push!(Error_k_B, abs(dd - did))
        push!(Error_k_C, abs(dd - did2))
        # push!(Error_k_C, ec)
        # push!(Error_k_D, ed)
    end
    # plot(k_array, log10.(Error_k_B))
    plot(k_array, log10.(Error_k_A))
    plot!(k_array, log10.(Error_k_B))
    plot!(k_array, log10.(Error_k_C))
end

#benchmark time
begin
    N_atoms = 100:100:10000
    time_cost = Vector{Float64}()
    for n_atoms in N_atoms
        q = 2 .* rand(n_atoms) .- 1.0
        q = q .- sum(q) / n_atoms

        x = 10.0 * rand(n_atoms)
        y = 10.0 * rand(n_atoms)
        z = 10.0 * rand(n_atoms)

        para = SoEwald2DPara((10.0, 10.0, 10.0), 1.0, 1.0, n_atoms);
        soepara = SoePara()
        iterpara = IterPara(n_atoms)
        adpara_dir = AdPara(Float64, n_atoms)
        adpara_soe = AdPara(Float64, n_atoms)
        K = (0.3, 0.4, 0.5)

        update_iterpara_z!(iterpara, z)
        # push!(time_cost, @belapsed energy_sum_k($K, $q, $x, $y, $z, $para, $soepara, $iterpara))
        push!(time_cost, @elapsed (for i in 1:10 energy_sum_k(K, q, x, y, z, para, soepara, iterpara) end))
    end
    plot(log10.(N_atoms), log10.(time_cost ./ 10), label = "time cost")
    plot!(xlabel = L"$\log(N)$", ylabel = L"$\log(T)$")
end