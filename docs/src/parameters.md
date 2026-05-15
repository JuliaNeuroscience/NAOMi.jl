```@meta
CurrentModule = NAOMi
```

# Parameters

Every pipeline stage is configured through a `Base.@kwdef mutable
struct`. Field names are preserved verbatim from the upstream MATLAB
`check_*_params.m` files (mixed `camelCase` / `snake_case`) so upstream
documentation maps one-to-one onto the Julia code. Defaults match
upstream.

Some structs carry *derived* fields whose defaults are sentinels
(`0`, `NaN`); [`finalize!`](@ref) resolves them in place.

## Volume and anatomy

```@docs
VolumeParams
NeuronParams
VasculatureParams
VasculatureNodeParams
DendriteParams
AxonParams
BackgroundParams
```

## Activity

```@docs
SpikeOptions
CalciumParams
```

## Optics and scanning

```@docs
PSFParams
PSFFastMask
TPMParams
ScanParams
NoiseParams
```

## Derived-field resolution

```@docs
finalize!
```
