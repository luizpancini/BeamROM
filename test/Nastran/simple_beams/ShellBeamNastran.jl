using BeamROM, LinearAlgebra
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

# 2D shell element static linear analysis with NASTRAN (SOL 101)

# Mesh ID (#1 = 50 x 10 elements
           #2 = 100 x 20 elements
           #3 = 200 x 40 elements)
mesh = 1

# Beam properties
L = 0.55
b = 60e-3
h = 2.25e-3
E = 71e9
ν = 0.33
G = E / (2 * (1 + ν))
ρ = 2.795e3

# Generalized displacements from the Nastran analysis (python scripts in /nastran)
T = nastran_static_response("shell", mesh; L=L, w=b, t=h, E=E, ν=ν, ρ=ρ)
# Displacements due to loads [m]
u_F1, u_F2, u_F3, u_M1, u_M2, u_M3 = (T[1:3, j] for j in 1:6)
# Rotation angles due to loads [rad]
θ_F1, θ_F2, θ_F3, θ_M1, θ_M2, θ_M3 = (T[4:6, j] for j in 1:6)

# Compute stiffness matrix
S = stiffness_matrix(L=L,u_F1=u_F1,u_F2=u_F2,u_F3=u_F3,u_M1=u_M1,u_M2=u_M2,u_M3=u_M3,θ_F1=θ_F1,θ_F2=θ_F2,θ_F3=θ_F3,θ_M1=θ_M1,θ_M2=θ_M2,θ_M3=θ_M3)

# Save the stiffness matrix next to the Nastran files of the mesh
write_stiffness_csv(joinpath(pkgdir(BeamROM), "nastran", "simple_beams", "shell_beam", "mesh$(mesh)", "shell_beam_stiffness_matrix.csv"), S)

# Comparison with expected stiffness matrix
# ------------------------------------------------------------------------------
# Derived geometric properties. Note: shear and torsional correction factors estimates based on http://dx.doi.org/10.17515/resm2015.19me0827 (or https://www.researchgate.net/publication/283665997_Shear_and_torsion_correction_factors_of_Timoshenko_beam_model_for_generic_cross_sections) for a rectangular cross-section with aspect ratio 60/2.25 ≈ 27
A = b*h
Iy = b*h^3/12
Iz = h*b^3/12
J = Iy+Iz
Kt = 6.5e-3
Ksy = Ksz = 5/6

# Expected stiffness matrix
S0 = diagm([E*A,G*A*Ksy,G*A*Ksz,G*J*Kt,E*Iy,E*Iz])

# Show results
println("Computed stiffness matrix:")
display(round.(S, digits=4))
println("Expected stiffness matrix:")
display(S0)
