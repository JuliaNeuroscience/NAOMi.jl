using NAOMi
using Random
using Statistics
using Test

@testset "spikes" begin

    @testset "sample_firing_rates :uniform returns identical rates" begin
        rng = MersenneTwister(1)
        so = SpikeOptions(K = 50, rate = 0.005, rate_dist = :uniform)
        r = sample_firing_rates(rng, so)
        @test length(r) == 50
        @test all(==(0.005), r)
    end

    @testset "sample_firing_rates :gamma respects clamp and mean" begin
        rng = MersenneTwister(2)
        base = 0.01
        so = SpikeOptions(K = 5000, rate = base, rate_dist = :gamma, alpha = 1.0)
        r = sample_firing_rates(rng, so)
        @test length(r) == 5000
        @test all(base/10 .≤ r .≤ base*10)
        # Empirical mean of exponential(scale=base), clamped to [base/10,10*base]
        # is within 25% of `base` (the clamp shaves a bit off the tail).
        @test mean(r) ≈ base rtol = 0.25
    end

    @testset "sample_firing_rates fallback symbol behaves like :uniform" begin
        rng = MersenneTwister(3)
        so = SpikeOptions(K = 10, rate = 0.02, rate_dist = :something_unknown)
        r = sample_firing_rates(rng, so)
        @test all(==(0.02), r)
    end

    @testset "generate_burst_spike_times returns correct shape & dtype" begin
        rng = MersenneTwister(4)
        so = SpikeOptions(K = 8, nt = 200, rate = 0.01,
                          rate_dist = :uniform, burst_mean = 0.0)
        S = generate_burst_spike_times(rng, so)
        @test size(S) == (8, 200)
        @test eltype(S) === Int
        @test all(x -> x == 0 || x == 1, S)
    end

    @testset "generate_burst_spike_times scales linearly with rate (no bursts)" begin
        rng = MersenneTwister(5)
        K, nt = 200, 2000
        rate = 0.002
        so = SpikeOptions(K = K, nt = nt, rate = rate,
                          rate_dist = :uniform, burst_mean = 0.0)
        S = generate_burst_spike_times(rng, so, fill(rate, K))
        # Average inter-spike interval = 1/rate = 500 samples → ~4 spikes/neuron
        # per nt=2000 samples. Per-neuron count ~ 4 ± few.
        avg_per_neuron = mean(sum(S, dims = 2))
        @test 2 < avg_per_neuron < 8

        # Doubling rate should roughly double the spike count.
        so2 = SpikeOptions(K = K, nt = nt, rate = 2*rate,
                           rate_dist = :uniform, burst_mean = 0.0)
        S2 = generate_burst_spike_times(MersenneTwister(5), so2, fill(2*rate, K))
        @test mean(sum(S2, dims = 2)) > 1.5 * avg_per_neuron
    end

    @testset "generate_burst_spike_times: bursts inflate spike count" begin
        rng = MersenneTwister(6)
        K, nt = 50, 5000
        rate = 0.001
        so0 = SpikeOptions(K = K, nt = nt, rate = rate,
                           rate_dist = :uniform, burst_mean = 0.0)
        sob = SpikeOptions(K = K, nt = nt, rate = rate,
                           rate_dist = :uniform, burst_mean = 4.0)
        S0 = generate_burst_spike_times(MersenneTwister(6), so0, fill(rate, K))
        Sb = generate_burst_spike_times(MersenneTwister(6), sob, fill(rate, K))
        @test sum(Sb) > sum(S0)
    end

    @testset "generate_burst_spike_times: reproducible with fixed RNG" begin
        so = SpikeOptions(K = 5, nt = 500, rate = 0.01, rate_dist = :uniform)
        S1 = generate_burst_spike_times(MersenneTwister(99), so)
        S2 = generate_burst_spike_times(MersenneTwister(99), so)
        @test S1 == S2
    end

    @testset "generate_burst_spike_times: rates length mismatch errors" begin
        so = SpikeOptions(K = 4)
        @test_throws ArgumentError generate_burst_spike_times(
            MersenneTwister(0), so, [1.0, 2.0])
    end

    @testset "bin_spike_trains: basic counts" begin
        # Events at t = 0.5, 1.5, 1.6 for neurons 1, 2, 2; dt = 1; T = 3
        evt = [0.5, 1.5, 1.6]
        evm = [1, 2, 2]
        S = bin_spike_trains(evt, evm, 3, 1.0, 3)
        @test size(S) == (3, 3)
        @test S[1, 1] == 1
        @test S[2, 2] == 2
        @test sum(S) == 3
    end

    @testset "bin_spike_trains: errors" begin
        @test_throws ArgumentError bin_spike_trains([1.0], [1, 2], 2, 1.0, 5)
        @test_throws ArgumentError bin_spike_trains([1.0], [99], 2, 1.0, 5)
        @test_throws ArgumentError bin_spike_trains([10.0], [1], 2, 1.0, 5)
    end

    @testset "bin_spike_trains: empty input" begin
        S = bin_spike_trains(Float64[], Int[], 3, 1.0, 4)
        @test size(S) == (3, 4)
        @test all(==(0), S)
    end

    @testset "sample_marked_point_process: homogeneous Poisson rate" begin
        # With constant CIF λ over [0, T], expect ≈ λ·T events.
        rng = MersenneTwister(11)
        λ, T = 5.0, 100.0
        cif(t, past_t, past_m) = λ
        times, marks = sample_marked_point_process(rng;
            cif = cif, timemax = T, nummax = 10_000)
        @test size(marks) == (length(times), 0)
        @test issorted(times)
        @test all(0 .≤ times .≤ T)
        @test length(times) ≈ λ * T rtol = 0.15
    end

    @testset "sample_marked_point_process: requires stopping criterion" begin
        @test_throws ArgumentError sample_marked_point_process(;
            cif = (t, pt, pm) -> 1.0)
    end

    @testset "sample_marked_point_process: nummax stops sampling" begin
        rng = MersenneTwister(12)
        times, _ = sample_marked_point_process(rng;
            cif = (t, pt, pm) -> 2.0,
            nummax = 7)
        @test length(times) == 7
    end

    @testset "sample_marked_point_process: marks are produced" begin
        rng = MersenneTwister(13)
        # Two-dim mark = (1.0, 2.0) for every event
        cif = (t, pt, pm) -> 1.0
        mkf = (t, pt, pm) -> (1.0, 2.0)
        times, marks = sample_marked_point_process(rng;
            cif = cif, mkf = mkf, markdim = 2, timemax = 5.0)
        @test size(marks, 2) == 2
        @test all(marks[:, 1] .== 1.0)
        @test all(marks[:, 2] .== 2.0)
    end

    @testset "sample_marked_point_process: monotone-decreasing CIF, custom cifmax" begin
        # Self-correcting process: rate is `λ0 · exp(-1·N)` where N is event count.
        # Decreases at events, constant between, so the default cifmax fallback
        # is correct.
        rng = MersenneTwister(14)
        λ0 = 10.0
        cif(t, past_t, past_m) = λ0 * exp(-length(past_t))
        times, _ = sample_marked_point_process(rng;
            cif = cif, timemax = 10.0, nummax = 1000)
        # Should eventually saturate; just a few events get through.
        @test length(times) < 100
        @test issorted(times)
    end
end
