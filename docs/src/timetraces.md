```@meta
CurrentModule = NAOMi
```

# Time traces

The time-trace stage turns simulated spikes into per-cell fluorescence
activity: spikes are drawn from a (possibly correlated) point process,
converted to calcium concentration through indicator dynamics, and
mapped to fluorescence. [`generate_time_traces`](@ref) is the top-level
entry point.

## Spike generation

```@docs
sample_firing_rates
generate_burst_spike_times
sample_marked_point_process
bin_spike_trains
```

## Calcium dynamics

```@docs
make_doub_exp_kernel
make_calcium_impulse
calcium_dynamics
fluorescence
```

## Top-level traces and correlation

```@docs
samp_small_world_mat
expression_variation
gen_correlated_spike_trains
generate_time_traces
```
