module VortexDetection3D

using VortexDistributions
using FourierGPE
using Plots
using LaTeXStrings
using JLD2
using FileIO
using Makie
using LinearAlgebra
using Interpolations

# density isosurface
export dense, densityfilm,

# functions in function.jld2
gradient_2D, gradient_3D, wps_v, wps_v_2D,
find_nearest, findidx_uniform, wps_array_plot,
gradient_3D_cent, grad_itp, wps_int, psi_itp,
vorticies_plot, find_closest_tuple,

# functions in algorithm
vortex_boundary_bottom, vortex_boundary_top,
vortex_link_bottom, vortex_link_top, sort_vorts_label,

# detection functions
find_vortices3D, find_vortices3D_v2, 

# psi Interpolations
psi_itp_2D

include("utils.jl")
include("detection_functions.jl")

@load joinpath(@__DIR__,"3dquenchslab_data.jld2") psi1 psi2 psi3
@load joinpath(@__DIR__,"3dquenchbox_data.jld2") psi_ring1 psi_ring2 psi_ring3
@load joinpath(@__DIR__,"3dquenchbox_data2.jld2") psi_box_25 psi_tubes
@load joinpath(@__DIR__,"3dquenchbox_data3.jld2") psi_ring4

export psi1, psi2, psi3,
psi_ring1, psi_ring2, psi_ring3,
psi_box_25, psi_tubes,
psi_ring4

end