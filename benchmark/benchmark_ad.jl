using BenchmarkTools, SoEwald2D, Plots

begin
    time_diff = []
    time_ad = []
    for n_atoms in 100:100:1000
        q = 2 .* rand(n_atoms) .- 1.0
        x = 10.0 * rand(n_atoms)
        y = 10.0 * rand(n_atoms)
        z = 10.0 * rand(n_atoms)

        U = [0.0]
        F_x = zeros(n_atoms)
        F_y = zeros(n_atoms)
        F_z = zeros(n_atoms)

        para = SoEwald2DPara((10.0, 10.0, 10.0), 1.0, 1.0, n_atoms);
        adpara = AdPara(Float64, n_atoms)

        t_d = @belapsed direct_sum!($q, $x, $y, $z, $para, $U), diff_direct_sum!($q, $x, $y, $z, $para, $F_x, $F_y, $F_z)
        t_ad = @belapsed ad_direct_sum!($q, $x, $y, $z, $para, $adpara)
        @show t_d, t_ad
        push!(time_diff, t_d)
        push!(time_ad, t_ad)
    end

    plot(n, time_diff, marker=:circle, label="manual diff")
    plot!(n, time_ad, marker=:circle, label="auto diff")
    xlabel!("n")
    ylabel!("time cost")
    savefig("ad_time.pdf")
end