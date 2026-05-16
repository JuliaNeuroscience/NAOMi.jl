```@meta
CurrentModule = NAOMi
```

# Volume

The volume stage samples a 3-D neural anatomy: blood vessels, cell
bodies, dendrites, axons, and neuropil background.
[`simulate_neural_volume`](@ref) runs the whole pipeline and returns a
[`NeuralVolume`](@ref).

## Top-level orchestration

```@docs
simulate_neural_volume
NeuralVolume
```

## Vasculature

```@docs
VesselNode
VesselEdge
gen_node
gen_conn
del_node!
nodes_to_conn
pseudo_rand_sample_2d
pseudo_rand_sample_3d
vessel_dijkstra
branch_grow_nodes!
grow_major_vessels!
grow_capillaries!
conn_to_vol!
simulate_blood_vessels
```

## Somata

```@docs
spiral_sample_sphere
teardrop_projection
generate_neural_body
sample_dense_neurons
generate_neural_volume
isolate_visible_somas
point_in_soma
masked_3d_gp
```

## Dendrites

```@docs
dendrite_dijkstra_grid
get_dendrite_path
dendrite_random_walk
dilate_dendrite_paths_all
grow_neuron_dendrites!
grow_apical_dendrites!
smooth_cell_body
set_cell_fluorescence
```

## Parallel dendrite growth

Dendrite growth dominates the cost of [`simulate_neural_volume`](@ref):
each neuron's basal and apical trees are traced with a grid Dijkstra
shortest-path search, and on a realistic volume that search runs for
tens of seconds per call. The searches for different neurons are almost
independent, so [`grow_neuron_dendrites!`](@ref) and
[`grow_apical_dendrites!`](@ref) can grow neurons in parallel — but not
quite for free.

Upstream NAOMi-Sim grows neurons **serially**, and each neuron's Dijkstra
search treats the dendrites of every previously-grown neuron as a soft
obstruction. Dendrites therefore tend to route *around* one another.
Growing all neurons in parallel breaks that coupling: every neuron's
search sees only the fixed somata, nuclei and vasculature, so dendrites
no longer avoid each other and overlap more.

The `couple_dendrites` keyword (accepted by [`simulate_neural_volume`](@ref)
and both dendrite-growth functions) selects between the two:

- `couple_dendrites = true` (the **default**) — upstream-faithful serial
  growth with inter-dendrite avoidance.
- `couple_dendrites = false` — parallel growth across
  `Threads.nthreads()` tasks. Several times faster (the dendrite stages
  of a 110 × 110 × 40 µm volume dropped from ~65 s to ~11 s in one
  benchmark), at the cost of the extra overlap quantified below.

Both modes are reproducible: each neuron draws from its own RNG seeded
deterministically from the caller's `rng`, and the two modes differ
*only* in the obstruction coupling.

### How much extra overlap?

Measured on a 110 × 110 × 40 µm volume with 49 neurons (identical
anatomy and seed), counting voxels occupied by the basal / local-apical
dendrites of two or more distinct neurons:

| Metric | `couple_dendrites = true` (default) | `couple_dendrites = false` |
|:-------|------------------------------------:|---------------------------:|
| Dendrite voxels | 290 278 | 252 426 |
| Voxels shared by ≥ 2 neurons | 2 258 (0.78 %) | 15 034 (5.96 %) |
| Most neurons through one voxel | 2 | 28 |

Decoupling raises the overlapping fraction roughly eight-fold (0.78 % →
5.96 %). It also lowers the total dendrite-voxel count: with coupling,
dendrites spread apart into more distinct voxels; without it they pile
onto the same low-cost corridors, so a single voxel can end up on the
path of dozens of neurons. The through-volume apical dendrites added by
[`grow_apical_dendrites!`](@ref) are far sparser (24 in this volume) and
showed no overlap in either mode — there the choice makes no difference.

Whether the extra overlap matters depends on the application. For
ground-truth movies where dendrites contribute diffuse neuropil
fluorescence it is usually negligible; if precise per-neuron dendrite
morphology matters, keep the default.

## Axons and neuropil background

```@docs
generate_axons
sort_axons
generate_bg_dendrites
```
