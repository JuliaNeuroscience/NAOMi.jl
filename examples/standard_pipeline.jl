# Reference end-to-end NAOMi pipeline.
#
# Julia translation of upstream `TPM_Simulation_Script_standard.m`
# (Copyright 2021 Alex Song, Adam Charles, MIT). Defaults to a small
# 30 × 30 × 20 µm volume so it runs in well under two minutes; the
# upstream script targets a 500 × 500 × 100 µm volume over 20 000 frames.
#
# Usage:
#   julia --project examples/standard_pipeline.jl [SX SY SZ [N_NEUR [NT]]]
#
#   SX SY SZ  volume size in microns          (default 30 30 20)
#   N_NEUR    number of neurons               (default 3 at the default
#             size; if a size is given without N_NEUR, the count is
#             derived from neuron density, so larger volumes are
#             automatically populated)
#   NT        number of imaging frames        (default 60)
#
# Larger volumes take proportionally longer (volume generation and
# scanning both scale with size). Output TIFFs land in a fresh temp
# directory whose path is printed at the end; set the NAOMI_OUTPUT_DIR
# environment variable to choose the location.

using NAOMi
using Random

# --- Command-line arguments ---------------------------------------------
# Returns (vol_sz::Vector{Int}, n_neur::Int, nt::Int). `n_neur == 0` is the
# "derive from density" sentinel resolved by `finalize!`.
function parse_pipeline_args(args)
    isempty(args) && return ([30, 30, 20], 3, 60)
    length(args) >= 3 ||
        error("Usage: standard_pipeline.jl [SX SY SZ [N_NEUR [NT]]]")
    nums = tryparse.(Int, args)
    any(isnothing, nums) &&
        error("All arguments must be integers; got $(args)")
    all(>(0), nums[1:3]) ||
        error("Volume size SX SY SZ must be positive; got $(nums[1:3])")
    vol_sz = nums[1:3]
    # n_neur == 0 is the "derive from density" sentinel resolved by `finalize!`.
    n_neur = length(nums) >= 4 ? nums[4] : 0
    n_neur >= 0 || error("N_NEUR must be ≥ 0 (0 derives from density); got $(n_neur)")
    nt = length(nums) >= 5 ? nums[5] : 60
    nt > 0 || error("NT must be positive; got $(nt)")
    return vol_sz, n_neur, nt
end

const SEED = 1
vol_sz, n_neur, nt = parse_pipeline_args(ARGS)

rng = MersenneTwister(SEED)

# Keep the volume centre comfortably below the brain surface.
vol_depth = max(15.0, 0.5 * vol_sz[3] + 5)
vol_params = VolumeParams(vol_sz=vol_sz, vol_depth=vol_depth, vres=2.0,
                          min_dist=12.0, N_neur=n_neur, N_bg=1)
finalize!(vol_params)                       # resolves N_neur if 0 was passed
vol_params.N_bg = 5 * vol_params.N_neur     # background processes scale with cells

neur_params = NeuronParams()
vasc_params = VasculatureParams(sourceFreq=200.0, vesFreq=[60.0, 80.0, 30.0])
dend_params = DendriteParams()
axon_params = AxonParams(maxlength=40.0, padsize=4, numbranches=1)
bg_params   = BackgroundParams(maxlength=40.0)
tpm_params  = TPMParams(pavg=20.0)
finalize!(tpm_params)
noise_params = NoiseParams()
scan_params  = ScanParams(scan_buff=4, motion=true)

println("Volume: $(vol_sz[1])×$(vol_sz[2])×$(vol_sz[3]) µm, " *
        "$(vol_params.N_neur) neurons, $(nt) frames")

# --- 1. Neural volume ----------------------------------------------------
# Grow dendrites in parallel when Julia was started with more than one
# thread (`julia -t auto …`); a single-threaded run keeps the slower,
# upstream-faithful serial growth. See the "Parallel dendrite growth"
# section of the docs for the speed/overlap tradeoff.
couple_dendrites = Threads.nthreads() == 1
println("Simulating neural volume (", couple_dendrites ? "serial" :
        "parallel, $(Threads.nthreads()) threads", " dendrites)...")
t0 = time()
neur_vol = simulate_neural_volume(vol_params, neur_params, vasc_params,
                                  dend_params, axon_params, bg_params;
                                  rng=rng, couple_dendrites)
println("  done in $(round(time() - t0; digits=2)) s")

# --- 2. Point-spread function -------------------------------------------
# Sampling is 1 / vres microns per voxel; a 15 × 15 × 11 voxel PSF.
psf, = gaussian_psf_na(tpm_params.nac, tpm_params.lambda,
                       1 / vol_params.vres, [15, 15, 11];
                       nidx=tpm_params.nidx)
psf = Float32.(psf ./ sum(psf))

# --- 3. Temporal activity -----------------------------------------------
println("Simulating temporal activity...")
spike_opts = SpikeOptions(K=vol_params.N_neur, nt=nt, dt=1 / 30,
                          rate=0.25, N_bg=axon_params.N_proc,
                          smod_flag=:independent)
traces = generate_time_traces(rng, spike_opts)
soma_act = Float32.(traces.soma)
dend_act = traces.dend === nothing ? copy(soma_act) : Float32.(traces.dend)
bg_act   = traces.bg   === nothing ?
    fill(0.5f0, axon_params.N_proc, nt) : Float32.(traces.bg)
neur_act = (soma=soma_act, dend=dend_act, bg=bg_act)

# --- 4. Scanning ---------------------------------------------------------
# Volume generation leaves a large amount of now-dead scratch on the heap
# (per-thread Dijkstra buffers). Collect it before scanning allocates the
# movie arrays, so the two stages' peaks do not stack up.
GC.gc()
println("Scanning volume ($(nt) frames)...")
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
# `cleanup=false` keeps the temp directory after the process exits — the
# default would delete the movies we just wrote.
out_dir = get(ENV, "NAOMI_OUTPUT_DIR",
              mktempdir(; prefix="naomi_demo_", cleanup=false))
isdir(out_dir) || mkpath(out_dir)
write_tpm_movie(joinpath(out_dir, "movie.tif"), mov)
write_tpm_movie(joinpath(out_dir, "movie_clean.tif"), mov_clean)

println()
println("Pipeline complete.")
println("  movie size      : $(size(mov))")
println("  ideal profiles  : $(size(comps))")
println("  ideal traces    : $(size(ideal_traces))")
println("  output directory: $out_dir")
