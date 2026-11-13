using BeamROM, LinearAlgebra, Printf, PrettyTables
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

# Global (entire wing) sectional stiffness matrix of Technion Pazy wing with 10 degrees of sweep (S10, Baseline version with no tip rod) from a static linear (SOL 101) analysis with NASTRAN. The wing is clamped at the root section and the six unit load cases are applied at the end of the reference axis, which is inclined by the sweep angle. The generalized displacements, the loads and the stiffness matrix are in the beam frame, which is rotated by the sweep angle about the normal axis: x1 along the reference axis (from the root to the end of the spar), x2 in the plane of the wing (towards the leading edge) and x3 normal to the wing

# Sweep angle [deg]: the sections (planes of constant y) of the Nastran model are shifted along x by y tan(sweep), so the reference axis is inclined by this angle from the y axis
sweep = 10

# Position of the reference axis, as a fraction of the chord (the same along the span, so that the axis is inclined by the sweep angle).
reference_axis = 0.447

# Generalized displacements from the Nastran static analysis of the Baseline without the tip rod (python scripts in /nastran), in the beam frame, and the length of the beam axis
T, L = nastran_pazy_swept_no_rod_static_response(sweep, "Baseline"; reference_axis)

# Displacements due to loads [m]
u_F1, u_F2, u_F3, u_M1, u_M2, u_M3 = (T[1:3, j] for j in 1:6)
# Rotation angles due to loads [rad]
θ_F1, θ_F2, θ_F3, θ_M1, θ_M2, θ_M3 = (T[4:6, j] for j in 1:6)

# Compute stiffness matrix
S = stiffness_matrix(L=L,u_F1=u_F1,u_F2=u_F2,u_F3=u_F3,u_M1=u_M1,u_M2=u_M2,u_M3=u_M3,θ_F1=θ_F1,θ_F2=θ_F2,θ_F3=θ_F3,θ_M1=θ_M1,θ_M2=θ_M2,θ_M3=θ_M3)

# Save the stiffness matrix next to the Nastran files (one file for each reference axis)
write_stiffness_csv(pazy_swept_no_rod_stiffness_file(sweep, "Baseline"; reference_axis), S)

# Show results
println("S10 Baseline without tip rod: length of the beam axis L = $L m")
println("Computed stiffness matrix:")
display(round.(S, digits=4))

# Modified stiffness matrix with updated shear stiffnesses
Smod = deepcopy(S)
Smod[3,3] = 1e6

# Generalized displacements from Nastran x from modified stiffness matrix
uθ = displacements_from_loads(loadCase="all",L=L,S=S)
uθMod = displacements_from_loads(loadCase="all",L=L,S=Smod)
println("Generalized displacements from original stiffness matrix:")
display(uθ)
println("Generalized displacements from modified stiffness matrix:")
display(uθMod)
println("Difference in generalized displacements due to modified shear stiffnesses:")
display(round.(uθMod .- uθ,digits=8))
println("Ratio of modified to original generalized displacements:")
display(round.(uθMod ./ uθ,digits=3))
