```@meta
CurrentModule = NAOMi
```

# Scanning

The scanning stage convolves the activity-modulated volume with the PSF
to produce a fluorescence movie, applies a Poisson–Gauss noise model and
motion, and extracts ground-truth "ideal" components for benchmarking
analysis pipelines. [`scan_volume`](@ref) is the top-level entry point.

## PSF FFT and single-frame scan

```@docs
psf_fft
single_scan
setup_scan_volume_frame
scan_volume_frame
ScanVolume
```

## Noise model

```@docs
poisson_gauss_noise
pixel_bleed
apply_noise_model
```

## Full scan and motion

```@docs
scan_volume
img_sub_row_shift
```

## Ideal components and ground truth

```@docs
calculate_ideal_comps
comps2ideals
times_from_profs
```
