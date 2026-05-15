using NAOMi
using Test

@testset "make_doub_exp_kernel" begin
    t_on, t_off, A, dt = 1.0, 5.0, 1.0, 0.01
    h = make_doub_exp_kernel(t_on, t_off, A, dt)
    @test h[1] == 0.0                              # kernel(0) = 0
    @test all(h .≥ -1e-12)
    # Peak should be near analytic loc_max = log((t_off+t_on)/t_off)/t_on
    loc_max = log((t_off + t_on) / t_off) / t_on
    peak_idx = argmax(h)
    @test isapprox((peak_idx - 1) * dt, loc_max; atol = 2 * dt)
    # Tail must drop below 1e-3 of peak by the last sample (definition of t_max)
    @test h[end] < 1e-3 * maximum(h)
    # Monotonically increasing before peak, decreasing after
    @test issorted(h[1:peak_idx])
    @test issorted(h[peak_idx:end]; rev = true)

    @test_throws ArgumentError make_doub_exp_kernel(-1.0, 1.0, 1.0, 0.01)
    @test_throws ArgumentError make_doub_exp_kernel(1.0, 1.0, 1.0, 0.0)
end

@testset "make_calcium_impulse" begin
    dt = 1 / 30
    h = make_calcium_impulse([0.5, 2.0]; dt = dt)
    @test length(h) == ceil(Int, 10 / dt)
    @test h[1] == 1.0
    # Two-pole closed form: h[n] = (p1^(n-1) − p2^(n-1)) / (p1 − p2) (shifted
    # to 1-based indexing).
    p1, p2 = exp(-0.5), exp(-2.0)
    expected = [(p1^n - p2^n) / (p1 - p2) for n in 1:length(h)]
    @test all(isapprox.(h, expected; atol = 1e-10))
    # Long-term decay should be dominated by the slower pole p1
    @test h[end] / h[end - 1] ≈ p1 atol = 1e-6

    @test_throws ArgumentError make_calcium_impulse([0.0, 0.0])
end

@testset "calcium_dynamics — sat_type :double" begin
    cp = CalciumParams(:GCaMP6f)        # sat_type default is :double
    @test cp.sat_type === :double
    K, nt = 4, 200
    S = zeros(K, nt)
    res = calcium_dynamics(S, cp)
    @test size(res.CB) == (K, nt)
    @test size(res.C) == (K, nt)
    @test size(res.F) == (K, nt)
    # No spikes → C stays at ca_rest (we initialised C[:,1] = ca_rest already)
    @test all(res.C[:, 1] .== cp.ca_rest)
    # Small slow drift comes from CB1/CB2 not starting at exact equilibrium
    # (upstream behaviour). Bound the drift loosely.
    @test maximum(abs.(res.C .- cp.ca_rest)) < 0.01 * cp.ca_rest
    # F at rest equals the Hill-equation baseline
    F_rest = res.F[1, end]
    @test isfinite(F_rest)
    @test F_rest > 1.0   # nonzero baseline contribution from ca_rest in the Hill curve

    # Impulse drives C up
    S2 = zeros(K, nt); S2[:, 10] .= 5e-5
    res2 = calcium_dynamics(S2, cp)
    @test maximum(res2.C[1, :]) > cp.ca_rest
    @test all(res2.F .≥ 1.0)
end

@testset "calcium_dynamics — sat_type :Ca_DE" begin
    cp = CalciumParams(:GCaMP6f)
    cp.sat_type = :Ca_DE
    K, nt = 2, 300
    S = zeros(K, nt); S[1, 20] = 1e-4
    res = calcium_dynamics(S, cp)
    @test size(res.C) == (K, nt)
    @test size(res.CB) == (K, nt)
    # Resting steady state when no spikes
    @test all(isapprox.(res.C[2, :], cp.ca_rest; atol = 1e-12))
    @test all(isapprox.(res.CB[2, :], cp.ca_rest; atol = 1e-12))
    # Spike drives the response above baseline, but eventually decays back
    @test maximum(res.CB[1, :]) > cp.ca_rest
    @test res.CB[1, end] ≈ cp.ca_rest atol = 1e-3 * cp.ca_amp
end

@testset "calcium_dynamics — sat_type :single" begin
    cp = CalciumParams(:GCaMP6f)
    cp.sat_type = :single
    K, nt = 2, 200
    S = zeros(K, nt); S[1, 5] = 1e-4
    res = calcium_dynamics(S, cp)
    @test size(res.C) == (K, nt)
    @test maximum(abs.(res.C[2, :] .- cp.ca_rest)) < 0.01 * cp.ca_rest
    @test maximum(res.C[1, :]) > cp.ca_rest
    @test all(isfinite, res.F)
end

@testset "calcium_dynamics — saturation cap" begin
    cp = CalciumParams(:GCaMP6f)
    cp.ca_sat = 0.5
    cap = cp.ca_dis * cp.ca_sat / (1 - cp.ca_sat)
    K, nt = 1, 100
    S = fill(0.0, K, nt); S[1, 5:10] .= 1.0       # huge drive
    res = calcium_dynamics(S, cp)
    @test maximum(res.C) ≤ cap + 1e-12
end

@testset "calcium_dynamics — over_samp" begin
    cp = CalciumParams(:GCaMP6f)
    K, nt = 2, 60
    S = zeros(K, nt); S[1, 10] = 1e-5
    r1 = calcium_dynamics(S, cp; over_samp = 1)
    r3 = calcium_dynamics(S, cp; over_samp = 3)
    @test size(r3.C) == (K, nt)
    @test size(r3.CB) == (K, nt)
end

@testset "calcium_dynamics — unknown sat_type" begin
    cp = CalciumParams(:GCaMP6f)
    cp.sat_type = :nonsense
    @test_throws ArgumentError calcium_dynamics(zeros(1, 3), cp)
end

@testset "fluorescence" begin
    CB = collect(range(1e-9, 5e-6; length = 50))
    for prot in (:GCaMP6, :GCaMP6f, :GCaMP6s, :GCaMP3, :OGB1, :OGB_1,
                 :GCaMP6_RS06, :GCaMP6rs06, :GCaMP6_RS09, :GCaMP6rs09,
                 :jGCaMP7f, :jGCaMP7s, :jGCaMP7b, :jGCaMP7c)
        F = fluorescence(CB, prot)
        @test all(isfinite, F)
        @test all(F .> 0)
        @test issorted(F)                # Hill curve is monotone in CB
    end
    # Unknown protein → warns and falls back to GCaMP6f values
    F_unknown = (@test_logs (:warn, r"Unknown protein") fluorescence(CB, :Bogus))
    F_known   = fluorescence(CB, :GCaMP6f)
    @test F_unknown == F_known
end
