```@meta
CurrentModule = NAOMi
```

# Optics

The optics stage builds the point-spread function used to scan the
volume: an analytic Gaussian PSF, optional Zernike-polynomial
aberrations applied at the back aperture, and Fresnel propagation of the
field through tissue.

## PSF kernels

```@docs
gaussian_psf
gaussian_psf_na
gaussian_beam_size
generate_gaussian_profile
```

## Zernike aberrations and back aperture

```@docs
zernike_polynomial
generate_zernike_weights
apply_zernike
generate_back_aperture
```

## Fresnel propagation

```@docs
fresnel_propagation_multi
group_z_project
width_estimate
width_estimate_3d
tpm_signal_scale
collection_mask
```
