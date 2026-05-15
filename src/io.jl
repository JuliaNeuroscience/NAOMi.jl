# TIFF movie I/O.
#
# Ported from upstream NAOMi-Sim files (Copyright 2021 Alex Song,
# Adam Charles, MIT):
#   - tifwrite.m / tiff_writer.m → write_tiff
#   - tifread.m  / tiff_reader.m → read_tiff
#   - tifwriteblock.m            → write_tiff_blocks
#   - write_TPM_movie.m          → write_tpm_movie
#
# Deferred (see ANALYSIS_PLAN.md, Chunk 19 inventory):
#   - tifinitialize.m / tifappend.m — frame-by-frame streaming append;
#     `write_tiff_blocks` covers the practical multi-file output need.
#   - make_avi.m — needs a plotting backend (`VideoIO`/`Makie`); the TIFF
#     path is the portable interchange format.
#   - saveSimulationParts.m — splits a MATLAB `.mat` into part files;
#     not meaningful for the Julia port (no monolithic `.mat`).

import TiffImages

export write_tiff, read_tiff, write_tiff_blocks, write_tpm_movie

# Block-file name: stem_00001.tif, stem_00002.tif, ...
_block_name(path::AbstractString, idx::Integer) =
    string(first(splitext(path)), "_", lpad(idx, 5, '0'), last(splitext(path)))

"""
    write_tiff(path::AbstractString, data::AbstractArray{<:Real};
               dtype::DataType=Float32)

Write a 2-D image or 3-D movie `data` to a (multi-page) TIFF at `path`.
Each `data[:, :, k]` becomes one page. Values are converted to `dtype`
(default `Float32`, matching upstream's `single` movies). Ports
`tifwrite.m` / `tiff_writer.m`.
"""
function write_tiff(path::AbstractString, data::AbstractArray{<:Real};
                    dtype::DataType=Float32)
    nd = ndims(data)
    (nd == 2 || nd == 3) ||
        throw(ArgumentError("image is not 2-D or 3-D (got $(nd)-D)"))
    arr = nd == 2 ? reshape(data, size(data, 1), size(data, 2), 1) : data
    gray = TiffImages.Gray.(dtype.(arr))
    TiffImages.save(path, gray)
    return path
end

"""
    read_tiff(path::AbstractString; frames=nothing) -> Array{Float32,3}

Read a TIFF stack at `path` into an `H × W × T` array of `Float32`.
`frames`, if given as a `(first, last)` tuple or range, truncates to that
inclusive page range. Ports `tifread.m` / `tiff_reader.m`.
"""
function read_tiff(path::AbstractString; frames=nothing)
    img = TiffImages.load(path)
    raw = Float32.(TiffImages.gray.(img))
    arr = ndims(raw) == 2 ?
        reshape(raw, size(raw, 1), size(raw, 2), 1) : Array(raw)
    isnothing(frames) && return arr
    lo, hi = first(frames), last(frames)
    (1 <= lo <= hi <= size(arr, 3)) ||
        throw(ArgumentError("frame range $(lo):$(hi) out of bounds " *
                            "for a $(size(arr, 3))-page stack"))
    return arr[:, :, lo:hi]
end

"""
    write_tiff_blocks(path::AbstractString, mov::AbstractArray{<:Real,3};
                      blocksize::Integer=2500, dtype::DataType=Float32)
        -> Vector{String}

Write a 3-D movie `mov` to a series of TIFF files, each holding at most
`blocksize` frames, named `stem_00001.tif`, `stem_00002.tif`, … (the
extension of `path` is reused). Returns the list of files written. Ports
`tifwriteblock.m`.
"""
function write_tiff_blocks(path::AbstractString, mov::AbstractArray{<:Real,3};
                           blocksize::Integer=2500, dtype::DataType=Float32)
    blocksize > 0 || throw(ArgumentError("blocksize must be positive"))
    T = size(mov, 3)
    n_blocks = cld(T, blocksize)
    written = String[]
    for b in 1:n_blocks
        lo = (b - 1) * blocksize + 1
        hi = min(b * blocksize, T)
        fname = _block_name(path, b)
        write_tiff(fname, @view(mov[:, :, lo:hi]); dtype=dtype)
        push!(written, fname)
    end
    return written
end

"""
    write_tpm_movie(path::AbstractString, mov::AbstractArray{<:Real,3};
                    blocksize::Integer=2500, dtype::DataType=Float32)
        -> Vector{String}

Write a simulated two-photon movie to disk, dispatching on the file
extension of `path`. Only `.tif` / `.tiff` is supported (via
[`write_tiff_blocks`](@ref)); upstream's `.fits` and `.mat` branches are
not ported. Ports `write_TPM_movie.m`.
"""
function write_tpm_movie(path::AbstractString, mov::AbstractArray{<:Real,3};
                         blocksize::Integer=2500, dtype::DataType=Float32)
    ext = lowercase(last(splitext(path)))
    if ext == ".tif" || ext == ".tiff"
        return write_tiff_blocks(path, mov; blocksize=blocksize, dtype=dtype)
    elseif ext == ".fits" || ext == ".mat"
        throw(ArgumentError("'$ext' output not ported; use a .tif path"))
    else
        throw(ArgumentError("unknown movie file type '$ext'; " *
                            "expected .tif"))
    end
end
