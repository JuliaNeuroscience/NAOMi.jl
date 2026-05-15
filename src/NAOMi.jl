module NAOMi

using Random
using Distributions

# Placeholder includes for the five-stage pipeline. Each is populated in a
# later chunk; see ANALYSIS_PLAN.md for the schedule.
include("params.jl")

include("timetraces/spikes.jl")
include("timetraces/calcium.jl")
include("timetraces/traces.jl")

include("optics/psf.jl")
include("optics/zernike.jl")
include("optics/propagation.jl")

include("volume/vasculature.jl")
include("volume/somata.jl")
include("volume/dendrites.jl")
include("volume/axons.jl")
include("volume/background.jl")
include("volume/volume.jl")

include("scanning/psf_fft.jl")
include("scanning/noise.jl")
include("scanning/scan.jl")
include("scanning/ideal.jl")

include("io.jl")

end # module NAOMi
