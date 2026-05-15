using NAOMi
using Test
using Random
using Statistics

@testset "vasculature helpers" begin
    @testset "VesselNode / VesselEdge defaults" begin
        n = NAOMi.gen_node()
        @test n.num == 0
        @test n.root == -1
        @test isempty(n.conn)
        @test n.type === :none
        n2 = NAOMi.gen_node(5, 2, [1, 3], [4.0, 5.0, 6.0], :surf, [1.0, 0.0])
        @test n2.num == 5
        @test n2.root == 2
        @test n2.conn == [1, 3]
        @test n2.pos == [4.0, 5.0, 6.0]
        @test n2.type === :surf
        @test n2.misc == [1.0, 0.0]

        e = NAOMi.gen_conn()
        @test e.start == 0
        @test isnan(e.weight)
        @test isempty(e.locs)
        @test e.misc === :none
        e2 = NAOMi.gen_conn(1, 2, 3.5, zeros(Int, 0, 3), :capp)
        @test e2.weight == 3.5
        @test e2.misc === :capp
    end

    @testset "del_node! clears slot and neighbours" begin
        nodes = [NAOMi.gen_node(1, 0, [2, 3], [0.0, 0.0, 0.0], :edge),
                 NAOMi.gen_node(2, 1, [1], [1.0, 0.0, 0.0], :surf),
                 NAOMi.gen_node(3, 1, [1], [0.0, 1.0, 0.0], :surf)]
        NAOMi.del_node!(nodes, 2)
        @test nodes[2].type === :none
        @test nodes[2].num == 0
        @test !(2 in nodes[1].conn)
        # 3's root pointer didn't reference 2, untouched
        @test nodes[3].root == 1
    end

    @testset "nodes_to_conn walks roots back from leaves" begin
        nodes = [NAOMi.gen_node(1, 0, [2], [0.0, 0.0, 0.0], :edge),
                 NAOMi.gen_node(2, 1, [1, 3], [0.0, 0.0, 1.0], :vert),
                 NAOMi.gen_node(3, 2, [2], [0.0, 0.0, 2.0], :capp, [4.0])]
        conn = NAOMi.nodes_to_conn(nodes)
        # one leaf (3), so a single chain of edges
        @test length(conn) == 2
        starts_ends = sort([(c.start, c.ends) for c in conn])
        @test starts_ends == [(2, 1), (3, 2)]
        # weight is leaf misc[1] = 4 carried back via sqrt-cumulation
        @test all(c -> isapprox(c.weight, 4.0; atol=1e-12), conn)
    end
end

@testset "pseudo_rand_sample" begin
    @testset "2-D positions stay in bounds and are distinct" begin
        rng = MersenneTwister(0)
        pos, _ = NAOMi.pseudo_rand_sample_2d((40, 40), 8;
                                             width=3.0, weight=0.5, rng=rng)
        @test size(pos, 2) == 2
        @test all(1 .<= pos .<= 40)
        # Mostly distinct samples thanks to Gaussian exclusion
        @test size(unique(pos; dims=1), 1) >= 6
    end

    @testset "3-D positions stay in bounds" begin
        rng = MersenneTwister(0)
        pos, _ = NAOMi.pseudo_rand_sample_3d((20, 20, 20), 6;
                                             width=2.0, weight=0.5, rng=rng)
        @test size(pos, 2) == 3
        @test all(1 .<= pos .<= 20)
    end

    @testset "near-corner sampling — bounds clipping" begin
        rng = MersenneTwister(0)
        # tight pdf forces samples; small grid forces near-edge sampling
        for _ in 1:5
            pos, _ = NAOMi.pseudo_rand_sample_2d((6, 6), 3; width=2.0, rng=rng)
            @test all(1 .<= pos .<= 6)
        end
    end
end

@testset "vessel_dijkstra" begin
    # Small reference graph (4 nodes, fully connected with known shortest paths).
    # Distances chosen so the unique min-cost paths are easy to verify.
    M = [0.0 1.0 4.0 Inf;
         1.0 0.0 2.0 5.0;
         4.0 2.0 0.0 1.0;
         Inf 5.0 1.0 0.0]
    d, pf = NAOMi.vessel_dijkstra(M, 1)
    @test d ≈ [0.0, 1.0, 3.0, 4.0]
    @test pf[1] == 0          # root has no predecessor
    @test pf[2] == 1
    @test pf[3] == 2
    @test pf[4] == 3
end

@testset "binary morphology" begin
    mask = falses(11, 11)
    mask[6, 6] = true
    NAOMi.dilate2d_disk!(mask, 2.0)
    # disk of radius 2 hits the four cardinal neighbours and the diagonals.
    @test mask[6, 6]
    @test mask[6, 8] && mask[8, 6] && mask[4, 6] && mask[6, 4]
    @test mask[5, 5]            # within √2 < 2
    @test !mask[6, 9]           # 3 > 2

    vol = falses(7, 7, 7)
    NAOMi.paint_ball3d!(vol, [4, 4, 4], 2.0, true)
    @test vol[4, 4, 4]
    @test vol[6, 4, 4] && vol[4, 6, 4] && vol[4, 4, 6]
    @test !vol[7, 7, 7]
end

@testset "simulate_blood_vessels — end-to-end" begin
    # Use a smallish volume; bump sourceFreq / vesFreq downward so that even
    # at 150 µm side length we get a few sources / surface / capp nodes —
    # otherwise rounding kills every count and the output is empty.
    vol = VolumeParams(vol_sz=[150, 150, 30], vol_depth=100.0, vres=2.0)
    finalize!(vol)
    vap = VasculatureParams(sourceFreq=400.0, vesFreq=[80.0, 100.0, 30.0])

    rng = MersenneTwister(0)
    neur_ves, _ = simulate_blood_vessels(vol, vap; rng=rng)
    @test ndims(neur_ves) == 3
    @test eltype(neur_ves) === Bool
    @test size(neur_ves) == (300, 300, 260)

    frac = sum(neur_ves) / length(neur_ves)
    @test 0.001 < frac < 0.5     # generous bracket; upstream ranges ~2-7%

    # Surface band (top ~30 µm = 60 voxels) should have non-zero density.
    surf = sum(@view neur_ves[:, :, 1:60])
    @test surf > 0
end

@testset "simulate_blood_vessels — reproducibility" begin
    vol = VolumeParams(vol_sz=[120, 120, 30], vol_depth=80.0, vres=2.0)
    finalize!(vol)
    vap = VasculatureParams(sourceFreq=400.0, vesFreq=[80.0, 100.0, 30.0])

    nv1, _ = simulate_blood_vessels(vol, vap; rng=MersenneTwister(7))
    nv2, _ = simulate_blood_vessels(vol, vap; rng=MersenneTwister(7))
    @test nv1 == nv2

    nv3, _ = simulate_blood_vessels(vol, vap; rng=MersenneTwister(8))
    @test nv1 != nv3
end
