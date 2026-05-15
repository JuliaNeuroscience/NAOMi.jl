```@meta
CurrentModule = NAOMi
```

# Getting started

This page walks through a complete simulation: build a neural volume,
synthesize a point-spread function, generate activity, scan the volume
into a movie, extract ground-truth components, and save the result to
TIFF. It mirrors `examples/standard_pipeline.jl`, scaled down to a small
30 × 30 × 20 µm volume so it runs in well under a minute.

## 1. Parameters

Every stage is configured through a `@kwdef` parameter struct (see
[Parameters](@ref)). Construct the ones the pipeline needs:

```julia
using NAOMi, Random

rng = MersenneTwister(1)

vol_params  = VolumeParams(vol_sz=[30, 30, 20], vol_depth=15.0, vres=2.0,
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
```

[`finalize!`](@ref) resolves derived fields (such as `N_neur` from a
density, or `phi` from the numerical aperture) that default to sentinel
values.

## 2. Neural volume

[`simulate_neural_volume`](@ref) runs the full anatomy pipeline —
vasculature, somata, dendrites, axons, neuropil — and returns a
[`NeuralVolume`](@ref):

```julia
neur_vol = simulate_neural_volume(vol_params, neur_params, vasc_params,
                                  dend_params, axon_params, bg_params;
                                  rng=rng)
```

## 3. Point-spread function

A Gaussian PSF parameterised by numerical aperture. The sampling is
`1 / vres` microns per voxel:

```julia
psf, = gaussian_psf_na(tpm_params.nac, tpm_params.lambda,
                       1 / vol_params.vres, [15, 15, 11];
                       nidx=tpm_params.nidx)
psf = Float32.(psf ./ sum(psf))
```

## 4. Temporal activity

[`generate_time_traces`](@ref) produces per-cell fluorescence traces
from simulated spikes and calcium dynamics:

```julia
spike_opts = SpikeOptions(K=vol_params.N_neur, nt=60, dt=1/30,
                          rate=0.25, N_bg=axon_params.N_proc,
                          smod_flag=:independent)
traces = generate_time_traces(rng, spike_opts)
neur_act = (soma = Float32.(traces.soma),
            dend = Float32.(traces.dend),
            bg   = Float32.(traces.bg))
```

## 5. Scanning

[`scan_volume`](@ref) convolves the activity-modulated volume with the
PSF frame by frame, applies motion and the Poisson–Gauss noise model,
and returns the movie. `return_clean=true` also yields the noise-free
movie:

```julia
mov, mov_clean = scan_volume(neur_vol, psf, neur_act, scan_params;
                             noise_params=noise_params,
                             tpm_params=tpm_params,
                             spike_opts=spike_opts, rng=rng,
                             return_clean=true)
```

## 6. Ground-truth components

[`calculate_ideal_comps`](@ref) produces per-cell spatial profiles, and
[`times_from_profs`](@ref) recovers per-cell traces from a movie given
those profiles:

```julia
comps, baseim, ideal =
    calculate_ideal_comps(neur_vol, psf, neur_act, scan_params;
                          noise_params=noise_params,
                          tpm_params=tpm_params,
                          spike_opts=spike_opts, rng=rng)
ideal_traces, = times_from_profs(mov_clean, comps; lambda=0)
```

## 7. Save to TIFF

[`write_tpm_movie`](@ref) writes the movie to disk:

```julia
write_tpm_movie("movie.tif", mov)
write_tpm_movie("movie_clean.tif", mov_clean)
```

The full script is at `examples/standard_pipeline.jl` and can be run
directly:

```
julia --project examples/standard_pipeline.jl
```
