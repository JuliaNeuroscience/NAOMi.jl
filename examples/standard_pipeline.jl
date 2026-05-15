# Reference end-to-end NAOMi pipeline.
#
# Julia translation of upstream `TPM_Simulation_Script_standard.m`
# (Copyright 2021 Alex Song, Adam Charles, MIT). Scaled down to a small
# 30 × 30 × 20 µm volume so it runs in well under two minutes; the
# upstream script targets a 500 × 500 × 100 µm volume over 20 000 frames.
#
# Run with:  julia --project examples/standard_pipeline.jl
# Output TIFFs land in a fresh temp directory (path printed at the end);
# set the NAOMI_OUTPUT_DIR environment variable to override.

using NAOMi
using Random

# --- Parameters ----------------------------------------------------------
const SEED  = 1
const NT    = 60          # number of imaging frames
const VOLSZ = [30, 30, 20]  # volume size in microns

rng = MersenneTwister(SEED)

vol_params  = VolumeParams(vol_sz=VOLSZ, vol_depth=15.0, vres=2.0,
                           min_dist=12.0, N_neur=3, N_bg=15)
finalize!(vol_params)
neur_params = NeuronParams()
vasc_params = VasculatureParams(sourceFreq=200.0, vesFreq=[60.0, 80.0, 30.0])
dend_params = DendriteParams()
axon_params = AxonParams(maxlength=40.0, padsize=4, numbranches=1)
bg_params   = BackgroundParams(maxlength=40.0)
tpm_params  = TPMParams(pavg=20.0)
finalize!(tpm_params)
noise_params = NoiseParams()
scan_params  = ScanParams(scan_buff=4, motion=true)

# --- 1. Neural volume ----------------------------------------------------
println("Simulating neural volume...")
t0 = time()
neur_vol = simulate_neural_volume(vol_params, neur_params, vasc_params,
                                  dend_params, axon_params, bg_params; rng=rng)
println("  done in $(round(time() - t0; digits=2)) s")

# --- 2. Point-spread function -------------------------------------------
# Sampling is 1 / vres microns per voxel; a 15 × 15 × 11 voxel PSF.
psf, = gaussian_psf_na(tpm_params.nac, tpm_params.lambda,
                       1 / vol_params.vres, [15, 15, 11];
                       nidx=tpm_params.nidx)
psf = Float32.(psf ./ sum(psf))

# --- 3. Temporal activity -----------------------------------------------
println("Simulating temporal activity...")
spike_opts = SpikeOptions(K=vol_params.N_neur, nt=NT, dt=1 / 30,
                          rate=0.25, N_bg=axon_params.N_proc,
                          smod_flag=:independent)
traces = generate_time_traces(rng, spike_opts)
soma_act = Float32.(traces.soma)
dend_act = traces.dend === nothing ? copy(soma_act) : Float32.(traces.dend)
bg_act   = traces.bg   === nothing ?
    fill(0.5f0, axon_params.N_proc, NT) : Float32.(traces.bg)
neur_act = (soma=soma_act, dend=dend_act, bg=bg_act)

# --- 4. Scanning ---------------------------------------------------------
println("Scanning volume ($(NT) frames)...")
t0 = time()
mov, mov_clean = scan_volume(neur_vol, psf, neur_act, scan_params;
                             noise_params=noise_params, tpm_params=tpm_params,
                             spike_opts=spike_opts, rng=rng,
                             return_clean=true)
println("  done in $(round(time() - t0; digits=2)) s")

# --- 5. Ideal ground-truth profiles + traces ----------------------------
println("Calculating ideal profiles...")
comps, baseim, ideal =
    calculate_ideal_comps(neur_vol, psf, neur_act, scan_params;
                          noise_params=noise_params, tpm_params=tpm_params,
                          spike_opts=spike_opts, rng=rng)
ideal_traces, = times_from_profs(mov_clean, comps; lambda=0)

# --- 6. Save outputs -----------------------------------------------------
out_dir = get(ENV, "NAOMI_OUTPUT_DIR", mktempdir(; prefix="naomi_demo_"))
isdir(out_dir) || mkpath(out_dir)
write_tpm_movie(joinpath(out_dir, "movie.tif"), mov)
write_tpm_movie(joinpath(out_dir, "movie_clean.tif"), mov_clean)

println()
println("Pipeline complete.")
println("  movie size      : $(size(mov))")
println("  ideal profiles  : $(size(comps))")
println("  ideal traces    : $(size(ideal_traces))")
println("  output directory: $out_dir")
