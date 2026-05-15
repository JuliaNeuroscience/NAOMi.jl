using NAOMi
using Test

@testset "params" begin

    @testset "VolumeParams defaults" begin
        vp = VolumeParams()
        @test vp.vol_sz == [100, 100, 50]
        @test vp.min_dist == 16.0
        @test vp.vres == 2.0
        @test vp.N_bg == 1_000_000
        @test vp.vol_depth == 200.0
        @test vp.verbose == 1
        @test vp.N_neur == 0
        @test vp.neur_density == 1e5
    end

    @testset "VolumeParams finalize! derives N_neur from density" begin
        vp = finalize!(VolumeParams())
        # default density 1e5 cells/mm³ on a 100x100x50 µm³ = 5e5 µm³ volume
        # → 1e5 * 5e5 / 1e9 = 50 cells
        @test vp.N_neur == 50
        @test vp.N_den ≈ 2e3 * (100 * 100) / 1e6
    end

    @testset "VolumeParams finalize! preserves user-supplied N_neur" begin
        vp = finalize!(VolumeParams(N_neur = 7))
        @test vp.N_neur == 7
        # neur_density should be back-computed
        @test vp.neur_density ≈ 1e9 * 7 / prod(vp.vol_sz)
    end

    @testset "VolumeParams finalize! rounds vol_sz[3] up to multiple of 10" begin
        vp = finalize!(VolumeParams(vol_sz = [100, 100, 33]))
        @test vp.vol_sz[3] == 40
    end

    @testset "VolumeParams keyword overrides" begin
        vp = VolumeParams(vol_sz = [50, 50, 20], vres = 4.0, verbose = 2)
        @test vp.vol_sz == [50, 50, 20]
        @test vp.vres == 4.0
        @test vp.verbose == 2
    end

    @testset "NeuronParams defaults match upstream" begin
        np = NeuronParams()
        @test np.n_samps == 200
        @test np.l_scale == 90.0
        @test np.p_scale == 1000.0
        @test np.avg_rad == 5.9
        @test np.nuc_rad == [5.65, 2.5]
        @test np.eccen == [0.35, 0.35, 0.5]
        @test np.neur_type === :pyr
    end

    @testset "VasculatureParams + nested node_params" begin
        vp = VasculatureParams()
        @test vp.flag == 1
        @test vp.ves_shift == [5.0, 15.0, 5.0]
        @test vp.vesSize == [15.0, 9.0, 2.0]
        @test vp.vesFreq == [125.0, 200.0, 50.0]
        @test vp.node_params isa VasculatureNodeParams
        @test vp.node_params.maxit == 25
        @test vp.node_params.lensc == 50.0
        @test vp.node_params.dirvar ≈ π / 8

        # Override the nested struct
        custom = VasculatureNodeParams(maxit = 100, mindist = 5.0)
        vp2 = VasculatureParams(node_params = custom)
        @test vp2.node_params.maxit == 100
        @test vp2.node_params.mindist == 5.0
    end

    @testset "DendriteParams defaults" begin
        dp = DendriteParams()
        @test dp.dtParams == [40.0, 150.0, 50.0, 1.0, 10.0]
        @test dp.atParams == [6.0, 5.0, 5.0, 5.0, 1.0]
        @test dp.rallexp == 1.5
        @test dp.dims == [60, 60, 60]
    end

    @testset "AxonParams defaults" begin
        ap = AxonParams()
        @test ap.flag == 1
        @test ap.fillweight == 100.0
        @test ap.maxel == 8
        @test ap.maxvoxel == 6
        @test ap.N_proc == 10
    end

    @testset "BackgroundParams defaults" begin
        bp = BackgroundParams()
        @test bp.maxel == 1
        @test bp.fillweight == 100.0
    end

    @testset "SpikeOptions defaults and symbols" begin
        so = SpikeOptions()
        @test so.K == 30
        @test so.dt ≈ 1/30
        @test so.nt == 1000
        @test so.dyn_type === :Ca_DE
        @test so.prot === :GCaMP6
        @test so.rate_dist === :gamma
        @test so.smod_flag === :hawkes
        @test so.burst_mean == 10.0
        @test so.spikeflag === true
    end

    @testset "CalciumParams default (GCaMP6)" begin
        cp = CalciumParams()
        @test cp.prot === :GCaMP6
        @test cp.ca_amp ≈ 76.1251
        @test cp.t_on ≈ 0.8535
        @test cp.t_off ≈ 98.6173
        @test cp.ext_rate ≈ 292.3
        @test cp.ca_dis ≈ 290e-9
    end

    @testset "CalciumParams protein-specific defaults" begin
        cp6f = CalciumParams(:GCaMP6f)
        @test cp6f.ca_amp ≈ 76.1251
        cp6s = CalciumParams(:GCaMP6s)
        @test cp6s.ca_amp ≈ 54.6943
        @test cp6s.t_on ≈ 0.4526
        @test cp6s.t_off ≈ 68.5461
        @test cp6s.ext_rate ≈ 299.0833
        cp3 = CalciumParams(:GCaMP3)
        @test cp3.ca_amp ≈ 0.05
        @test cp3.t_on == 1.0
        cp7 = CalciumParams(:GCaMP7)
        @test cp7.ca_amp ≈ 230.917
        @test cp7.t_on ≈ 0.020137
    end

    @testset "CalciumParams keyword overrides protein default" begin
        cp = CalciumParams(:GCaMP6s, ca_amp = 99.0)
        @test cp.prot === :GCaMP6s
        @test cp.ca_amp == 99.0
        # Untouched kinetic constants still get the GCaMP6s default
        @test cp.t_off ≈ 68.5461
    end

    @testset "CalciumParams unknown protein falls back to GCaMP6 defaults" begin
        cp = CalciumParams(:JuliaFP)
        @test cp.ca_amp ≈ 76.1251
        @test cp.t_on ≈ 0.8535
    end

    @testset "NoiseParams defaults match upstream effective values" begin
        np = NoiseParams()
        @test np.mu == 100.0
        @test np.sigma == 2300.0
        @test np.sigma0 == 2.7
        @test np.darkcount == 0.05
        # Effective upstream values (the second duplicate block is dead code)
        @test np.bleedp == 0.3
        @test np.bleedw == 0.4
    end

    @testset "ScanParams defaults" begin
        sp = ScanParams()
        @test sp.scan_buff == 10
        @test sp.motion === true
        @test sp.scan_avg == 2
        @test sp.sfrac == 2
        @test sp.verbose == 1
    end

    @testset "PSFParams defaults and FastMask" begin
        psf = PSFParams()
        @test psf.NA == 0.6
        @test psf.objNA == 0.8
        @test psf.n == 1.35
        @test psf.lambda ≈ 0.92
        @test psf.psf_sz == [20.0, 20.0, 50.0]
        @test psf.scatter_sz == [0.51, 1.56, 4.52, 14.78]
        @test psf.type === :gaussian
        @test psf.scaling === :two_photon
        @test psf.hemoabs ≈ 0.00674 * log(10)
        @test psf.fastmask === true
        @test psf.FM isa PSFFastMask
        @test psf.FM.sampling == 10.0
        @test psf.FM.fineSamp == 2.0
        @test psf.FM.ss == 1.0
    end

    @testset "TPMParams defaults and phi derivation" begin
        tp = TPMParams()
        @test tp.nidx == 1.33
        @test tp.nac == 0.8
        @test isnan(tp.phi)
        @test tp.eta == 0.6
        @test tp.pavg == 40.0
        @test tp.lambda ≈ 0.92

        finalize!(tp)
        @test !isnan(tp.phi)
        # Upstream formula: 0.8 * ((1 - sqrt(1 - (nac/nidx)^2)) / 2) * 0.4
        expected = 0.8 * ((1 - sqrt(1 - (0.8/1.33)^2)) / 2) * 0.4
        @test tp.phi ≈ expected
    end

    @testset "TPMParams user-supplied phi survives finalize!" begin
        tp = finalize!(TPMParams(phi = 0.123))
        @test tp.phi == 0.123
    end

    @testset "finalize! is idempotent" begin
        vp = VolumeParams()
        finalize!(vp)
        n1 = vp.N_neur
        finalize!(vp)
        @test vp.N_neur == n1
    end

    @testset "every parameter struct has a default constructor" begin
        # If any field is missing a default, this fails immediately on `T()`.
        for T in (VolumeParams, NeuronParams, VasculatureParams,
                  VasculatureNodeParams, DendriteParams, AxonParams,
                  BackgroundParams, SpikeOptions, CalciumParams,
                  NoiseParams, PSFParams, PSFFastMask, ScanParams,
                  TPMParams)
            @test T() isa T
        end
    end
end
