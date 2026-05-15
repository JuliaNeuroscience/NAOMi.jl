using NAOMi
using Test

@testset "write_tiff / read_tiff — 3-D round trip" begin
    mov = rand(Float32, 14, 11, 6)
    path = tempname() * ".tif"
    write_tiff(path, mov)
    back = read_tiff(path)
    @test size(back) == size(mov)
    @test eltype(back) === Float32
    @test back ≈ mov
end

@testset "write_tiff / read_tiff — 2-D round trip" begin
    img = rand(Float32, 9, 13)
    path = tempname() * ".tif"
    write_tiff(path, img)
    back = read_tiff(path)
    @test size(back) == (9, 13, 1)
    @test back[:, :, 1] ≈ img
end

@testset "read_tiff — frame range truncation" begin
    mov = reshape(Float32.(1:(8 * 8 * 10)), 8, 8, 10)
    path = tempname() * ".tif"
    write_tiff(path, mov)
    sub = read_tiff(path; frames=(3, 6))
    @test size(sub) == (8, 8, 4)
    @test sub ≈ mov[:, :, 3:6]
    @test_throws "out of bounds" read_tiff(path; frames=(5, 99))
end

@testset "write_tiff — rejects non-image arrays" begin
    @test_throws "not 2-D or 3-D" write_tiff(tempname() * ".tif",
                                             rand(Float32, 3, 3, 3, 3))
end

@testset "write_tiff_blocks — splits a movie into block files" begin
    mov = rand(Float32, 10, 12, 7)
    path = tempname() * ".tif"
    files = write_tiff_blocks(path, mov; blocksize=3)
    @test length(files) == 3                       # 7 frames → 3 + 3 + 1
    @test all(isfile, files)
    # Reassembling the block files recovers the movie.
    recon = cat((read_tiff(f) for f in files)...; dims=3)
    @test size(recon) == size(mov)
    @test recon ≈ mov
    @test_throws "blocksize must be positive" write_tiff_blocks(path, mov;
                                                                blocksize=0)
end

@testset "write_tpm_movie — extension dispatch" begin
    mov = rand(Float32, 6, 6, 4)
    files = write_tpm_movie(tempname() * ".tif", mov; blocksize=10)
    @test length(files) == 1
    @test read_tiff(files[1]) ≈ mov
    @test_throws "not ported" write_tpm_movie(tempname() * ".mat", mov)
    @test_throws "unknown movie file type" write_tpm_movie(tempname() * ".xyz",
                                                           mov)
end

@testset "write_tiff — dtype conversion to Float64" begin
    mov = rand(Float32, 5, 5, 3)
    path = tempname() * ".tif"
    write_tiff(path, mov; dtype=Float64)
    back = read_tiff(path)
    @test size(back) == size(mov)
    @test back ≈ mov
end
