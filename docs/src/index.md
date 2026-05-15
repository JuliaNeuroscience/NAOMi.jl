```@meta
CurrentModule = NAOMi
```

# NAOMi.jl

NAOMi.jl is a Julia port of the MATLAB **NAOMi-Sim** simulator for
two-photon calcium imaging — a synthetic ground-truth generator for
evaluating calcium-imaging analysis methods. It samples a 3-D neural
anatomy (somata, dendrites, axons, vasculature, neuropil), simulates
calcium dynamics, propagates light through tissue, and scans the
resulting volume to produce realistic noisy fluorescence movies.

The pipeline is organised into five stages, each with its own API
section:

| Stage | What it does | API page |
|:------|:-------------|:---------|
| Time traces | Spikes → calcium → fluorescence activity | [Time traces](@ref) |
| Optics | Gaussian PSF, Zernike aberrations, Fresnel propagation | [Optics](@ref) |
| Volume | Vasculature, somata, dendrites, axons, neuropil | [Volume](@ref) |
| Scanning | PSF convolution, noise, motion, ideal components | [Scanning](@ref) |
| I/O | TIFF movie read/write | [I/O](@ref) |

All five stages share a set of `@kwdef` parameter structs documented on
the [Parameters](@ref) page.

See [Getting started](@ref) for an end-to-end example.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/JuliaNeuroscience/NAOMi.jl")
```

## Citation

If you use NAOMi.jl in published work, please cite the original NAOMi
paper:

> Song, A., Gauthier, J. L., Pillow, J. W., Tank, D. W., and Charles,
> A. S. (2021). Neural anatomy and optical microscopy (NAOMi)
> simulation for evaluating calcium imaging methods. *Journal of
> Neuroscience Methods*, 358, 109173.
> [doi:10.1016/j.jneumeth.2021.109173](https://doi.org/10.1016/j.jneumeth.2021.109173)

The upstream MATLAB implementation by Alex Song and Adam Charles is at
<https://bitbucket.org/adamshch/naomi_sim>. Attribution and licensing
details are recorded in `NOTICE.md`.

## Index

```@index
```
