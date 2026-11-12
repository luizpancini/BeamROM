using BeamROM, LinearAlgebra, Printf, PrettyTables
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

# Global (entire wing) sectional stiffness matrix of Technion Pazy wing without the tip rod from a static linear (SOL 101) analysis with NASTRAN. The wing is clamped at the root section and the six unit load cases are applied at the end of the reference axis.

# Position of the reference axis, as a fraction of the chord
reference_axis = 0.4653

# Generalized displacements from the Nastran static analysis (python scripts in /nastran), in the beam frame (x1 = y along the span from the root to the end of the spar, x2 = -x, x3 = z), and the length of the beam
T, L = nastran_pazy_no_rod_static_response(; reference_axis)
# Displacements due to loads [m]
u_F1, u_F2, u_F3, u_M1, u_M2, u_M3 = (T[1:3, j] for j in 1:6)
# Rotation angles due to loads [rad]
θ_F1, θ_F2, θ_F3, θ_M1, θ_M2, θ_M3 = (T[4:6, j] for j in 1:6)

# Compute stiffness matrix
S = stiffness_matrix(L=L,u_F1=u_F1,u_F2=u_F2,u_F3=u_F3,u_M1=u_M1,u_M2=u_M2,u_M3=u_M3,θ_F1=θ_F1,θ_F2=θ_F2,θ_F3=θ_F3,θ_M1=θ_M1,θ_M2=θ_M2,θ_M3=θ_M3)

# Save the stiffness matrix next to the Nastran files (one file for each reference axis)
write_stiffness_csv(pazy_no_rod_stiffness_file(; reference_axis), S)

# Comparison with the stiffness matrix of the spar alone (60 x 2.25 mm aluminum plate, see OneDBeamNastran.jl for the correction factors), as a reference. The wing also has skin, ribs and the tip, so the stiffnesses are expected to be close, but not equal
b = 60e-3
h = 2.25e-3
E = 71e9
ν = 0.33
G = E / (2 * (1 + ν))
ρ = 2.795e3
A = b*h
Iy = b*h^3/12
Iz = h*b^3/12
J = Iy+Iz
Kt = 6.5e-3
Ksy = Ksz = 5/6
S0 = diagm([E*A,G*A*Ksy,G*A*Ksz,G*J*Kt,E*Iy,E*Iz])

# Show results
println("Computed stiffness matrix:")
display(round.(S, digits=4))
println("Spar stiffness matrix:")
display(S0)

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