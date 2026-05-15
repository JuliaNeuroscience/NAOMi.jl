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

## Axons and neuropil background

```@docs
generate_axons
sort_axons
generate_bg_dendrites
```
