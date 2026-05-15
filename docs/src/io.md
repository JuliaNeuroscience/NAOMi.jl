```@meta
CurrentModule = NAOMi
```

# I/O

Simulated movies are written to and read from multi-page TIFF stacks via
[`TiffImages.jl`](https://github.com/tlnagy/TiffImages.jl).
[`write_tpm_movie`](@ref) is the recommended entry point for saving a
finished movie.

```@docs
write_tiff
read_tiff
write_tiff_blocks
write_tpm_movie
```
