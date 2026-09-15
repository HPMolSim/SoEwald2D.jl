@testset "framework-free plan API (Task 3)" begin
    # Every bound below is derived from a measured deviation on the same
    # configuration, with the margin stated. No bound is a guess, and none is
    # so loose that the assertion cannot fail.

    reldev(a, b) = (d = max(abs(a), abs(b)); d > 0 ? abs(a - b) / d : 0.0)

    # ---- reference configuration ------------------------------------------
    n = 20
    L = (15.0, 15.0, 15.0)
    ϵ_0 = 1.0
    α = 0.3
    s = 2.0
    r_c = s / α         # r_c = 6.6667 < min(Lx, Ly) / 2 = 7.5
    k_c = 2 * s * α
    @test r_c < min(L[1], L[2]) / 2

    Random.seed!(20260918)
    poses = [SVector(rand() * L[1], rand() * L[2], 1.0 + rand() * (L[3] - 2.0)) for _ in 1:n]
    charges = [(isodd(i) ? 1.0 : -1.0) * (1 + 0.02i) for i in 1:n]
    charges .-= sum(charges) / n          # exactly neutral

    long_plan = SoEwald2DLongPlan(ϵ_0, L, s, α, n, k_c, SoePara(); parallel = false)
    short_plan = SoEwald2DShortPlan(ϵ_0, L, s, α, n, r_c)

    @testset "long-range energy and force against the direct sums" begin
        El = SoEwald2D.energy(long_plan, poses, charges)
        # measured rel dev 3.56e-14 against direct_sum -> 1e-11 is a 280x margin
        @test isapprox(El, direct_sum(long_plan, poses, charges), rtol = 1e-11)
        # measured abs dev 2.50e-16 between the erfc and SOE direct sums
        # -> 1e-13 is a 400x margin
        @test isapprox(direct_sum(long_plan, poses, charges),
                       real(soe_direct_sum(long_plan, poses, charges, long_plan.soepara)),
                       atol = 1e-13)

        Fl = SoEwald2D.force(long_plan, poses, charges)
        ref = diff_direct_sum(long_plan, poses, charges)
        # measured: max rel dev 3.33e-13, max abs dev 1.86e-15, |F|max 9.14e-3
        # -> rtol 1e-10 is a 300x margin, atol 1e-13 a 70x margin
        for i in 1:n, d in 1:3
            @test isapprox(Fl[i][d], ref[i][d], rtol = 1e-10, atol = 1e-13)
        end
    end

    @testset "force! and force agree, and force! fills rather than accumulates" begin
        for plan in (short_plan, long_plan)
            F = SoEwald2D.force(plan, poses, charges)
            buf = [SVector(1.0, 2.0, 3.0) for _ in 1:n]   # deliberately dirty
            SoEwald2D.force!(buf, plan, poses, charges)
            # exact: same code path, same order, and force! must overwrite the
            # junk above rather than add to it
            @test buf == F
            @test SoEwald2D.force!(buf, plan, poses, charges) === buf
        end
    end

    @testset "queries never mutate poses or charges" begin
        poses_0 = deepcopy(poses)
        charges_0 = deepcopy(charges)
        SoEwald2D.energy(short_plan, poses, charges)
        SoEwald2D.energy(long_plan, poses, charges)
        SoEwald2D.force(short_plan, poses, charges)
        SoEwald2D.force(long_plan, poses, charges)
        direct_sum(long_plan, poses, charges)
        soe_direct_sum(long_plan, poses, charges, long_plan.soepara)
        diff_direct_sum(long_plan, poses, charges)
        @test poses == poses_0
        @test charges == charges_0
    end

    @testset "a supplied neighbor_list is candidate-only: its reported r is ignored" begin
        # Phase 1 of this project produced a wrong-SIGNED Ewald2D energy
        # (+0.0238 against a true -0.1539) by trusting a neighbour list's own
        # reported distance, which was in-plane for a 2-D finder. So the list
        # here carries a nonsense distance (-1.0) in the third slot, in the
        # same (i, j) order the O(n^2) loop uses; the results must be
        # BIT-identical, not merely close.
        full = [(i, j, -1.0) for i in 1:n for j in (i + 1):n]
        @test SoEwald2D.energy(short_plan, poses, charges; neighbor_list = full) ==
              SoEwald2D.energy(short_plan, poses, charges)
        @test SoEwald2D.force(short_plan, poses, charges; neighbor_list = full) ==
              SoEwald2D.force(short_plan, poses, charges)
    end

    @testset "_min_image_slab returns the FULL 3-D distance, not the in-plane one" begin
        # The single most likely way to break this package: copying
        # QuasiEwald's `_min_image_q2d`, which returns the in-plane rho_sq
        # because quasi-2D Ewald needs that. SoEwald2D's kernel is
        # erfc(alpha*r)/r in the TRUE separation.
        #
        # Two particles in one column, 2.0 apart in z and identical in x and y.
        # In-plane: rho_sq = 0, which the coincident-pair guard would then skip,
        # leaving only the two self terms. In 3-D: r = 2.0 < r_c, a real pair
        # contribution. So this testset fails loudly for an in-plane helper.
        L2 = (20.0, 20.0, 20.0)
        α2 = 0.2
        s2 = 1.0
        r_c2 = s2 / α2      # r_c = 5.0 < min(Lx, Ly) / 2 = 10.0
        @test r_c2 < min(L2[1], L2[2]) / 2
        sp2 = SoEwald2DShortPlan(1.0, L2, s2, α2, 2, r_c2)
        col = [SVector(4.0, 5.0, 6.0), SVector(4.0, 5.0, 8.0)]
        q2 = [1.0, -1.0]

        _, _, r_sq = SoEwald2D._min_image_slab(col[1], col[2], L2)
        @test r_sq == 4.0                        # dz^2, not the in-plane 0

        self_only = sum(SoEwald2D.SoEwald2D_Es_self(q2[i], α2) for i in 1:2) / (4π * 1.0)
        pair_term = q2[1] * q2[2] * erfc(α2 * 2.0) / 2.0 / (4π * 1.0)
        E = SoEwald2D.energy(sp2, col, q2)
        # measured |E - self - pair| = 3.47e-18 -> atol 1e-16 is a 29x margin
        @test isapprox(E, self_only + pair_term, atol = 1e-16)
        # and the pair term is not a rounding-level addition, so an in-plane
        # helper (which would drop it) cannot pass the line above by accident
        @test abs(pair_term) > 1e-3

        # q1 = +1 at z = 6 attracted to q2 = -1 at z = 8, so F_1z > 0, and the
        # pair force is exactly antisymmetric with no in-plane component.
        F = SoEwald2D.force(sp2, col, q2)
        @test F[1][3] > 0
        @test F[1][3] == -F[2][3]
        @test F[1][1] == 0.0 && F[1][2] == 0.0
    end

    @testset "force is -dE/dr: every component, sign included" begin
        # A sign error is invisible to any magnitude-only test, and the
        # pre-decoupling long-range path returned the energy GRADIENT
        # (`acceleration -= Point(Fx, Fy, Fz) / mass`), so the flip to a force
        # is exactly the kind of change that needs pinning down. `isapprox`
        # against the central difference asserts magnitude AND sign at once;
        # the explicit product test below makes the sign intent unmissable.

        @testset "long plan, (s, α) = ($sl, $al)" for (sl, al) in ((1.0, 0.5), (2.0, 0.4), (3.0, 0.3))
            nl = 6
            Ll = (10.0, 10.0, 10.0)
            Random.seed!(1)
            pl = [SVector(rand() * Ll[1], rand() * Ll[2], 1.0 + rand() * (Ll[3] - 2.0)) for _ in 1:nl]
            ql = [(isodd(i) ? 1.0 : -1.0) * (1 + 0.02i) for i in 1:nl]
            ql .-= sum(ql) / nl
            plan = SoEwald2DLongPlan(1.0, Ll, sl, al, nl, 2 * sl * al, SoePara(); parallel = false)
            F = SoEwald2D.force(plan, pl, ql)
            h = 1e-4
            # measured max rel dev over the three (s, α) pairs at h = 1e-4:
            # 3.61e-9 (max abs dev 5.85e-12) -> rtol 1e-6 is a 280x margin.
            # h was swept over 1e-4 ... 1e-7: no component disagreed at any h,
            # so no two-dimensional accuracy sweep was needed.
            for i in 1:nl, d in 1:3
                e = SVector{3, Float64}(ntuple(k -> k == d ? h : 0.0, 3))
                pp = copy(pl); pp[i] = pl[i] + e
                pm = copy(pl); pm[i] = pl[i] - e
                fd = -(SoEwald2D.energy(plan, pp, ql) - SoEwald2D.energy(plan, pm, ql)) / (2h)
                @test isapprox(F[i][d], fd, rtol = 1e-6, atol = 1e-12)
            end
            # the sign, stated on its own for the largest component
            i_max, d_max = argmax([abs(F[i][d]) for i in 1:nl, d in 1:3]).I
            e = SVector{3, Float64}(ntuple(k -> k == d_max ? h : 0.0, 3))
            pp = copy(pl); pp[i_max] = pl[i_max] + e
            pm = copy(pl); pm[i_max] = pl[i_max] - e
            dE = (SoEwald2D.energy(plan, pp, ql) - SoEwald2D.energy(plan, pm, ql)) / (2h)
            @test F[i_max][d_max] * dE < 0        # force opposes the gradient
        end

        @testset "short plan" begin
            # A compact cluster, so that every pair sits at ~0.8 << r_c = 5.0
            # and no finite difference straddles the cutoff, where
            # erfc(α r_c)/r_c is discontinuous.
            Ls = (20.0, 20.0, 20.0)
            αs = 0.2
            ss = 1.0
            r_cs = ss / αs      # r_c = 5.0 < min(Lx, Ly) / 2 = 10.0
            @test r_cs < min(Ls[1], Ls[2]) / 2
            Random.seed!(7)
            ps = [SVector(10.0 + 0.8 * (rand() - 0.5), 10.0 + 0.8 * (rand() - 0.5),
                          10.0 + 0.8 * (rand() - 0.5)) for _ in 1:6]
            qs = [(isodd(i) ? 1.0 : -1.0) * (1 + 0.03i) for i in 1:6]
            @test maximum(sqrt(sum(abs2, ps[i] - ps[j])) for i in 1:6, j in 1:6 if i < j) < 0.5 * r_cs
            plan = SoEwald2DShortPlan(1.0, Ls, ss, αs, 6, r_cs)
            F = SoEwald2D.force(plan, ps, qs)
            h = 1e-5
            # measured max rel dev 6.11e-9 at h = 1e-5 (|F|max 3.67)
            # -> rtol 1e-6 is a 160x margin
            for i in 1:6, d in 1:3
                e = SVector{3, Float64}(ntuple(k -> k == d ? h : 0.0, 3))
                pp = copy(ps); pp[i] = ps[i] + e
                pm = copy(ps); pm[i] = ps[i] - e
                fd = -(SoEwald2D.energy(plan, pp, qs) - SoEwald2D.energy(plan, pm, qs)) / (2h)
                @test isapprox(F[i][d], fd, rtol = 1e-6, atol = 1e-12)
            end
            i_max, d_max = argmax([abs(F[i][d]) for i in 1:6, d in 1:3]).I
            e = SVector{3, Float64}(ntuple(k -> k == d_max ? h : 0.0, 3))
            pp = copy(ps); pp[i_max] = ps[i_max] + e
            pm = copy(ps); pm[i_max] = ps[i_max] - e
            dE = (SoEwald2D.energy(plan, pp, qs) - SoEwald2D.energy(plan, pm, qs)) / (2h)
            @test F[i_max][d_max] * dE < 0
        end
    end

    @testset "r = 0 is finite and contributes exactly zero" begin
        # The `iszero(r_sq)` half of the short-range guard. Replacing
        # `position_check3D`'s all-zero sentinel with a plain `r_sq >= r_c^2`
        # test would drop it, and `SoEwald2D_Fs_pair` divides by the pair
        # separation -- so a coincident pair becomes 0/0 = NaN and poisons
        # every other particle's force through the pairwise accumulation. This
        # exact substitution shipped that NaN in QuasiEwald.
        L0 = (20.0, 20.0, 20.0)
        α0 = 0.2
        s0 = 1.0
        r_c0 = s0 / α0      # r_c = 5.0 < min(Lx, Ly) / 2 = 10.0
        @test r_c0 < min(L0[1], L0[2]) / 2
        q0 = [1.0, -0.5, 0.75]

        # (a) two particles at literally the same point;
        # (b) a pair whose x separation is an exact multiple of Lx at equal
        #     y and z, so the minimum image wraps to exactly zero -- which any
        #     lattice initialisation using inclusive endpoints will produce.
        cases = Dict(
            "coincident" => [SVector(4.0, 5.0, 6.0), SVector(4.0, 5.0, 6.0), SVector(6.0, 5.0, 7.0)],
            "wrapped to zero" => [SVector(0.0, 5.0, 6.0), SVector(20.0, 5.0, 6.0), SVector(2.0, 5.0, 7.0)],
        )
        for (label, p0) in cases
            @testset "$label" begin
                _, _, r01 = SoEwald2D._min_image_slab(p0[1], p0[2], L0)
                @test r01 == 0.0            # the guard really is reached
                plan = SoEwald2DShortPlan(1.0, L0, s0, α0, 3, r_c0)

                E = SoEwald2D.energy(plan, p0, q0)
                F = SoEwald2D.force(plan, p0, q0)
                @test isfinite(E)
                @test all(isfinite, (F[i][d] for i in 1:3, d in 1:3))

                # "contributes exactly zero", not merely "is finite": the
                # total must equal the same sum with the (1,2) pair term
                # contributing nothing. Accumulated in the query's own order --
                # pairs (1,2), (1,3), (2,3), then the three self terms -- so
                # the comparison can be exact rather than approximate; a
                # different summation order costs 2 ULPs here and would force a
                # tolerance that a genuinely wrong pair term could hide behind.
                E_ref = zero(Float64)
                E_ref += SoEwald2D._short_pair_energy(plan, p0, q0, 1, 2, r_c0^2)
                E_ref += SoEwald2D._short_pair_energy(plan, p0, q0, 1, 3, r_c0^2)
                E_ref += SoEwald2D._short_pair_energy(plan, p0, q0, 2, 3, r_c0^2)
                for i in 1:3
                    E_ref += SoEwald2D.SoEwald2D_Es_self(q0[i], α0)
                end
                @test SoEwald2D._short_pair_energy(plan, p0, q0, 1, 2, r_c0^2) == 0.0
                @test E == E_ref / (4π * 1.0)

                f13 = SoEwald2D._short_pair_force(plan, p0, q0, 1, 3, r_c0^2) / (4π * 1.0)
                f23 = SoEwald2D._short_pair_force(plan, p0, q0, 2, 3, r_c0^2) / (4π * 1.0)
                @test F[1] == f13
                @test F[2] == f23
                @test F[3] == -f13 - f23
                # and the pair that WAS kept is not itself zero, so the three
                # equalities above are not vacuously true
                @test any(!iszero, f13)
            end
        end

        # And one coincident pair must not poison the rest of the array: a
        # fourth, well-separated particle keeps exactly the force it would have
        # had on its own.
        p0 = [SVector(4.0, 5.0, 6.0), SVector(4.0, 5.0, 6.0), SVector(6.0, 5.0, 7.0), SVector(9.0, 5.0, 6.5)]
        q4 = [1.0, -0.5, 0.75, -1.25]
        plan4 = SoEwald2DShortPlan(1.0, L0, s0, α0, 4, r_c0)
        F4 = SoEwald2D.force(plan4, p0, q4)
        @test all(isfinite, (F4[i][d] for i in 1:4, d in 1:3))
        @test all(!isnan, (F4[i][d] for i in 1:4, d in 1:3))
    end

    @testset "SoEwald2DShortPlan rejects r_c >= min(Lx, Ly) / 2" begin
        # The single nearest in-plane image is the whole story only below half
        # the box; at or above it a second image is also inside the cutoff and
        # is silently dropped. Nothing downstream can detect that, so the
        # constructor is the only place it can be caught for a standalone
        # caller.
        @test_throws ArgumentError SoEwald2DShortPlan(1.0, (10.0, 10.0, 10.0), 1.0, 0.2, 4, 5.0)   # r_c == min/2
        @test_throws ArgumentError SoEwald2DShortPlan(1.0, (10.0, 12.0, 10.0), 1.0, 0.2, 4, 5.5)   # r_c > min/2 = 5.0
        @test SoEwald2DShortPlan(1.0, (10.0, 10.0, 10.0), 1.0, 0.25, 4, 4.0).r_c == 4.0            # 4.0 < 5.0
        msg = try
            SoEwald2DShortPlan(1.0, (10.0, 10.0, 10.0), 1.0, 0.2, 4, 5.0)
            ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("r_c < min(Lx, Ly) / 2", msg)
    end

    @testset "plans accept Point and NTuple positions with no conversion layer" begin
        # The interchange type is SVector{3,T}, but the kernels index only
        # p[1]/p[2]/p[3], so ExTinyMD's Point and a plain NTuple must give the
        # identical answer. Bit-identical, not approximately: there is no
        # conversion step for a rounding difference to hide in.
        as_point = [Point(p[1], p[2], p[3]) for p in poses]
        as_tuple = [(p[1], p[2], p[3]) for p in poses]
        for alt in (as_point, as_tuple)
            @test SoEwald2D.energy(short_plan, alt, charges) == SoEwald2D.energy(short_plan, poses, charges)
            @test SoEwald2D.energy(long_plan, alt, charges) == SoEwald2D.energy(long_plan, poses, charges)
            @test SoEwald2D.force(short_plan, alt, charges) == SoEwald2D.force(short_plan, poses, charges)
            @test SoEwald2D.force(long_plan, alt, charges) == SoEwald2D.force(long_plan, poses, charges)
        end
    end

    @testset "the rbm (random batch) path runs, and its serial branch is alive" begin
        # Each plan carries its own `rng = MersenneTwister(123)` and every rbm
        # query draws from it, so a plan that has already been queried is NOT
        # comparable with a fresh one. Every comparison below therefore gets
        # its own freshly constructed plan.
        mk(par) = SoEwald2DLongPlan(ϵ_0, L, s, α, n, k_c, SoePara(); rbm = true, rbm_p = 12, parallel = par)

        E1 = SoEwald2D.energy(mk(false), poses, charges)
        @test isfinite(E1)
        @test E1 == SoEwald2D.energy(mk(true), poses, charges)

        F1 = SoEwald2D.force(mk(false), poses, charges)
        @test all(isfinite, (F1[i][d] for i in 1:n, d in 1:3))
        # `parallel = false` with `rbm = true` used to be dead code: `for i in
        # indice` referred to a plan FIELD that was never bound as a local, so
        # the branch raised an unconditional UndefVarError. It is `1:rbm_p`
        # now, which is what the parallel branch uses and what `random_k` is
        # indexed by, so the two branches must agree exactly. If the serial
        # branch were still broken this line would error rather than fail.
        @test SoEwald2D.force(mk(true), poses, charges) == F1
    end
end
