using NAOMi
using Test
using Random
using Statistics

@testset "poisson_gauss_noise — empirical mean/variance match theory" begin
    nps = NoiseParams()
    rng = MersenneTwister(0)
    # Sweep across several clean-intensity levels and check matched
    # mean and variance against the analytic predictions.
    for lam_in in (1.0f0, 10.0f0, 50.0f0)
        clean = fill(lam_in, 200, 200)
        noisy = poisson_gauss_noise(clean, nps; rng=rng)
        emp_mean = mean(noisy)
        emp_var  = var(noisy)
        lam = lam_in + Float32(nps.darkcount)
        th_mean = lam * nps.mu
        th_var  = lam * nps.sigma + lam * nps.mu^2 + nps.sigma0^2
        @test isapprox(emp_mean, th_mean; rtol=0.05)
        @test isapprox(emp_var,  th_var;  rtol=0.10)
    end
end

@testset "poisson_gauss_noise — zero clean gives near-zero output" begin
    nps = NoiseParams(darkcount=0.0)
    rng = MersenneTwister(0)
    clean = zeros(Float32, 100, 100)
    noisy = poisson_gauss_noise(clean, nps; rng=rng)
    # Only electronic noise contributes; mean ≈ mu0, std ≈ sigma0 (rounded).
    @test abs(mean(noisy)) <= 1.0
    @test std(noisy) <= 1.5 * nps.sigma0
end

@testset "pixel_bleed — preserves total signal" begin
    rng = MersenneTwister(0)
    frame = rand(rng, Float32, 30, 30) .* 100
    bled = pixel_bleed(frame, 0.3, 0.4; rng=MersenneTwister(0))
    # Each pixel donates x_bleed·frame to its successor; net total
    # changes only at the very first pixel (which has no predecessor)
    # and may also drift slightly from float-precision rounding.
    @test isapprox(sum(bled), sum(frame); rtol=0.05)
end

@testset "pixel_bleed — p == 0 is a no-op" begin
    frame = rand(Float32, 10, 10)
    @test pixel_bleed(frame, 0.0, 0.4) ≈ Float32.(frame)
end

@testset "apply_noise_model — movie shape preserved" begin
    nps = NoiseParams()
    rng = MersenneTwister(0)
    mov = fill(5.0f0, 20, 20, 4)
    out = apply_noise_model(mov, nps; rng=rng)
    @test size(out) == size(mov)
    @test eltype(out) === Float32
    # Per-frame mean is positive (signal preserved).
    for k in 1:4
        @test mean(@view(out[:, :, k])) > 0
    end
end
