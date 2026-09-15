@testset "core works without ExTinyMD" begin
    # The requirement this whole phase exists for: SoEwald2D's core query API
    # must work with ExTinyMD never loaded. A @testset inside the normal suite
    # does not prove that on its own -- test/runtests.jl itself does
    # `using ExTinyMD` (for the adapter tests), so by the time this testset
    # runs ExTinyMD is already loaded in *this* process. The only way to prove
    # the core does not need it is to run in a fresh process that never
    # imports it.
    script = """
    using SoEwald2D, StaticArrays, SpecialFunctions
    const EXTINYMD = Base.PkgId(Base.UUID("fec76197-d59f-46dd-a0ed-76a83c21f7aa"), "ExTinyMD")
    @assert !haskey(Base.loaded_modules, EXTINYMD)

    n = 24
    L = (12.0, 12.0, 8.0)
    poses = [SVector(rand()*L[1], rand()*L[2], 1.0 + rand()*(L[3]-2.0)) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    ϵ_0 = 1.0
    α = 0.4
    s = 2.0
    r_c = s / α      # r_c = 5.0 < min(Lx,Ly)/2 = 6.0
    k_c = 2 * s * α

    sp = SoEwald2DShortPlan(ϵ_0, L, s, α, n, r_c)
    lp = SoEwald2DLongPlan(ϵ_0, L, s, α, n, k_c, SoePara(); parallel = false)

    Es = SoEwald2D.energy(sp, poses, charges)
    El = SoEwald2D.energy(lp, poses, charges)
    Fs = SoEwald2D.force(sp, poses, charges)
    Fl = SoEwald2D.force(lp, poses, charges)
    @assert isfinite(Es) && isfinite(El)
    @assert all(isfinite, Fs[i][d] for i in 1:n, d in 1:3)
    @assert all(isfinite, Fl[i][d] for i in 1:n, d in 1:3)

    # force! into a caller-owned buffer, and the NTuple position form
    buf = [SVector(0.0, 0.0, 0.0) for _ in 1:n]
    @assert SoEwald2D.force!(buf, sp, poses, charges) === buf
    @assert buf == Fs
    @assert SoEwald2D.energy(sp, [(p[1], p[2], p[3]) for p in poses], charges) == Es

    # the accuracy references are framework-free too
    @assert isfinite(direct_sum(lp, poses, charges))
    @assert isfinite(real(soe_direct_sum(lp, poses, charges, lp.soepara)))
    dd = diff_direct_sum(lp, poses, charges)
    @assert all(isfinite, dd[i][d] for i in 1:n, d in 1:3)

    # ... and still not loaded, after all of that
    @assert !haskey(Base.loaded_modules, EXTINYMD)
    print("OK")
    """
    out = read(`$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $script`, String)
    @test out == "OK"

    # Verify the check above is not vacuous: it must actually fail if ExTinyMD
    # is loaded first.
    #
    # `@test !success(proc)` alone is NOT enough, and that is the whole point:
    # it is true for any non-zero exit -- ExTinyMD missing from the test
    # environment, an unsatisfiable resolve, a failed precompile, a typo in the
    # heredoc. Every one of those would make this check pass while proving
    # nothing about whether the `@assert` on `Base.loaded_modules` is
    # load-bearing. So the subprocess's streams are captured and the failure is
    # pinned to the intended cause: a marker printed after `using` and before
    # the assert must appear on stdout (so the loads all succeeded), the marker
    # after the assert must NOT appear (so the assert is what stopped it), and
    # stderr must carry the AssertionError naming the guard's own expression.
    script_contaminated = """
    using ExTinyMD, SoEwald2D, StaticArrays
    print("LOADS_OK;")
    @assert !haskey(Base.loaded_modules, Base.PkgId(
        Base.UUID("fec76197-d59f-46dd-a0ed-76a83c21f7aa"), "ExTinyMD"))
    print("ASSERT_DID_NOT_TRIP;")
    """
    outfile, errfile = tempname(), tempname()
    proc = run(pipeline(
        `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $script_contaminated`;
        stdout = outfile, stderr = errfile), wait = false)
    wait(proc)
    out_text = read(outfile, String)
    err_text = read(errfile, String)
    rm(outfile, force = true)
    rm(errfile, force = true)

    @test !success(proc)
    # `using ExTinyMD, SoEwald2D, StaticArrays` really did succeed, so the
    # non-zero exit is not a missing package, a resolve failure or a
    # precompilation error.
    @test occursin("LOADS_OK;", out_text)
    # and execution stopped at the assert, not after it.
    @test !occursin("ASSERT_DID_NOT_TRIP;", out_text)
    # and it stopped for exactly the intended reason.
    @test occursin("AssertionError", err_text)
    @test occursin("loaded_modules", err_text)

    # The dispatcher's two error cases. `Base.get_extension` returns `nothing`
    # both when ExTinyMD was never loaded and when it IS loaded but the
    # extension failed to precompile; reporting the second as the first tells
    # the user to run a `using` they have already run, while the real error has
    # scrolled past. Only the first case is reachable from a test (the second
    # needs a deliberately broken extension), so it is the one asserted -- with
    # the "loaded but the extension failed" wording explicitly excluded, which
    # is what pins the branch.
    script_no_extinymd = """
    using SoEwald2D
    for f in (SoEwald2DShortInteraction, SoEwald2DLongInteraction)
        try
            f(1.0, (10.0, 10.0, 10.0), 1.0, 0.25, 2, 4.0)
            print("NO_ERROR;")
        catch e
            print(sprint(showerror, e), ";")
        end
    end
    """
    msg = read(`$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $script_no_extinymd`, String)
    @test count("requires ExTinyMD to be loaded", msg) == 2   # both names, not just one
    @test occursin("SoEwald2DExTinyMDExt", msg)
    @test occursin("SoEwald2DShortPlan", msg)                 # says what to use instead
    @test !occursin("NO_ERROR", msg)
    @test !occursin("retry_load_extensions", msg)             # that is the *other* branch
end
