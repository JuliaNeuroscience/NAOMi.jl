# NAOMi

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaNeuroscience.github.io/NAOMi.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaNeuroscience.github.io/NAOMi.jl/dev/)
[![Build Status](https://github.com/JuliaNeuroscience/NAOMi.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaNeuroscience/NAOMi.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaNeuroscience/NAOMi.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaNeuroscience/NAOMi.jl)

NAOMi.jl is a Julia port (in progress) of the MATLAB **NAOMi-Sim** simulator
for two-photon calcium imaging — a synthetic ground-truth generator that
samples a 3-D neural anatomy (somata, dendrites, axons, vasculature,
neuropil), simulates calcium dynamics, propagates light through tissue, and
scans the resulting volume to produce realistic noisy fluorescence movies.

The port is being implemented in chunks; see [`ANALYSIS_PLAN.md`](ANALYSIS_PLAN.md)
for the schedule and current status. Attribution and licensing for the
upstream work are recorded in [`NOTICE.md`](NOTICE.md).

If you use NAOMi.jl in published work, please cite the original NAOMi paper:

> Song, A., Gauthier, J. L., Pillow, J. W., Tank, D. W., and Charles, A. S.
> (2021). Neural anatomy and optical microscopy (NAOMi) simulation for
> evaluating calcium imaging methods. *Journal of Neuroscience Methods*,
> 358, 109173.

Upstream MATLAB code: <https://bitbucket.org/adamshch/naomi_sim>.
