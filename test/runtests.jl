using NAOMi
using Test

@testset "NAOMi.jl" begin
    include("test_params.jl")
    include("timetraces/test_spikes.jl")
    include("timetraces/test_calcium.jl")
    include("timetraces/test_traces.jl")
end
