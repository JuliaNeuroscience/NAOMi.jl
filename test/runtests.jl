using NAOMi
using Test

@testset "NAOMi.jl" begin
    include("test_params.jl")
    include("timetraces/test_spikes.jl")
    include("timetraces/test_calcium.jl")
    include("timetraces/test_traces.jl")
    include("optics/test_psf.jl")
    include("optics/test_zernike.jl")
    include("optics/test_propagation.jl")
    include("volume/test_vasculature.jl")
    include("volume/test_somata.jl")
end
