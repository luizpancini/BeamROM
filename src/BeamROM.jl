module BeamROM

using ForwardDiff
using LinearAlgebra
using Plots
using Printf
using Revise

include("StiffnessCore.jl")
include("StiffnessUtilities.jl")
include("InertiaCore.jl")
include("InertiaUtilities.jl")

export stiffness_matrix, stiffness_matrix_derivatives, stiffness_matrix_normalized_derivatives, displacements_from_loads, discretized_stiffness_matrices, clamped_segment_stiffness_matrices, least_squares_stiffness_matrices, isotropic_stiffness_matrix, anisotropic_stiffness_matrix, write_stiffness_csv, write_discretized_stiffness_csv

export beam_mass_model, straight_axis, read_gfem, lump_gfem, lump_by_grid, condense_to_beam, rigid_body_mass_matrix, write_beam_mass_model, write_beam_mass_csv, pazy_beam_mass_model, spar_reference_axis, rib_stations

end