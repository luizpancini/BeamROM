using BeamROM, LinearAlgebra

# 1D beam element static linear analysis with Ansys for a steel beam with rectangular cross-sectional dimensions 50 x 2 mm

# Beam length [m]
L = 1

# Displacements due to loads [m]
u_F1 = [5e-8; 0; 0]
u_F2 = [0; 8.015437197e-5; 0]
u_F3 = [0; 0; 5.000014603e-2]
u_M1 = [0; 0; 0]
u_M2 = [0; 0; -7.499998808e-2]
u_M3 = [0; 1.19999997e-4; 0]

# Rodrigues parameters (may be taken as approximately equal to the angles of rotation themselves, if small) due to loads
θ_F1 = [0; 0; 0]
θ_F2 = [0; 0; 1.2e-4]
θ_F3 = [0; -7.499998549e-2; 0]
θ_M1 = [9.798876027e-2; 0; 0]
θ_M2 = [0; 1.499999772e-1; 0]
θ_M3 = [0; 0; 2.4e-4]

# Compute stiffness matrix
S = stiffness_matrix(L=L,u_F1=u_F1,u_F2=u_F2,u_F3=u_F3,u_M1=u_M1,u_M2=u_M2,u_M3=u_M3,θ_F1=θ_F1,θ_F2=θ_F2,θ_F3=θ_F3,θ_M1=θ_M1,θ_M2=θ_M2,θ_M3=θ_M3)

# Show result
println("Computed stiffness matrix:")
display(S)

# Stiffness matrix sensitivity to displacements
J,S11_uF1,S11_θF1,S11_uF2,S11_θF2,S11_uF3,S11_θF3,S11_uM1,S11_θM1,S11_uM2,S11_θM2,S11_uM3,S11_θM3,S12_uF1,S12_θF1,S12_uF2,S12_θF2,S12_uF3,S12_θF3,S12_uM1,S12_θM1,S12_uM2,S12_θM2,S12_uM3,S12_θM3,S13_uF1,S13_θF1,S13_uF2,S13_θF2,S13_uF3,S13_θF3,S13_uM1,S13_θM1,S13_uM2,S13_θM2,S13_uM3,S13_θM3,S14_uF1,S14_θF1,S14_uF2,S14_θF2,S14_uF3,S14_θF3,S14_uM1,S14_θM1,S14_uM2,S14_θM2,S14_uM3,S14_θM3,S15_uF1,S15_θF1,S15_uF2,S15_θF2,S15_uF3,S15_θF3,S15_uM1,S15_θM1,S15_uM2,S15_θM2,S15_uM3,S15_θM3,S16_uF1,S16_θF1,S16_uF2,S16_θF2,S16_uF3,S16_θF3,S16_uM1,S16_θM1,S16_uM2,S16_θM2,S16_uM3,S16_θM3,S22_uF1,S22_θF1,S22_uF2,S22_θF2,S22_uF3,S22_θF3,S22_uM1,S22_θM1,S22_uM2,S22_θM2,S22_uM3,S22_θM3,S23_uF1,S23_θF1,S23_uF2,S23_θF2,S23_uF3,S23_θF3,S23_uM1,S23_θM1,S23_uM2,S23_θM2,S23_uM3,S23_θM3,S24_uF1,S24_θF1,S24_uF2,S24_θF2,S24_uF3,S24_θF3,S24_uM1,S24_θM1,S24_uM2,S24_θM2,S24_uM3,S24_θM3,S25_uF1,S25_θF1,S25_uF2,S25_θF2,S25_uF3,S25_θF3,S25_uM1,S25_θM1,S25_uM2,S25_θM2,S25_uM3,S25_θM3,S26_uF1,S26_θF1,S26_uF2,S26_θF2,S26_uF3,S26_θF3,S26_uM1,S26_θM1,S26_uM2,S26_θM2,S26_uM3,S26_θM3,S33_uF1,S33_θF1,S33_uF2,S33_θF2,S33_uF3,S33_θF3,S33_uM1,S33_θM1,S33_uM2,S33_θM2,S33_uM3,S33_θM3,S34_uF1,S34_θF1,S34_uF2,S34_θF2,S34_uF3,S34_θF3,S34_uM1,S34_θM1,S34_uM2,S34_θM2,S34_uM3,S34_θM3,S35_uF1,S35_θF1,S35_uF2,S35_θF2,S35_uF3,S35_θF3,S35_uM1,S35_θM1,S35_uM2,S35_θM2,S35_uM3,S35_θM3,S36_uF1,S36_θF1,S36_uF2,S36_θF2,S36_uF3,S36_θF3,S36_uM1,S36_θM1,S36_uM2,S36_θM2,S36_uM3,S36_θM3,S44_uF1,S44_θF1,S44_uF2,S44_θF2,S44_uF3,S44_θF3,S44_uM1,S44_θM1,S44_uM2,S44_θM2,S44_uM3,S44_θM3,S45_uF1,S45_θF1,S45_uF2,S45_θF2,S45_uF3,S45_θF3,S45_uM1,S45_θM1,S45_uM2,S45_θM2,S45_uM3,S45_θM3,S46_uF1,S46_θF1,S46_uF2,S46_θF2,S46_uF3,S46_θF3,S46_uM1,S46_θM1,S46_uM2,S46_θM2,S46_uM3,S46_θM3,S55_uF1,S55_θF1,S55_uF2,S55_θF2,S55_uF3,S55_θF3,S55_uM1,S55_θM1,S55_uM2,S55_θM2,S55_uM3,S55_θM3,S56_uF1,S56_θF1,S56_uF2,S56_θF2,S56_uF3,S56_θF3,S56_uM1,S56_θM1,S56_uM2,S56_θM2,S56_uM3,S56_θM3,S66_uF1,S66_θF1,S66_uF2,S66_θF2,S66_uF3,S66_θF3,S66_uM1,S66_θM1,S66_uM2,S66_θM2,S66_uM3,S66_θM3 = stiffness_matrix_derivatives(L=L,u_F1=u_F1,u_F2=u_F2,u_F3=u_F3,u_M1=u_M1,u_M2=u_M2,u_M3=u_M3,θ_F1=θ_F1,θ_F2=θ_F2,θ_F3=θ_F3,θ_M1=θ_M1,θ_M2=θ_M2,θ_M3=θ_M3)

# Comparison with expected stiffness matrix
# ------------------------------------------------------------------------------
# Beam properties
b = 50e-3
h = 2e-3
E = 200e9
G = 76.92e9

# Derived geometric properties. Note: shear and torsional correction factors estimates based on http://dx.doi.org/10.17515/resm2015.19me0827 (or https://www.researchgate.net/publication/283665997_Shear_and_torsion_correction_factors_of_Timoshenko_beam_model_for_generic_cross_sections) for a rectangular cross-section with aspect ratio 50/2 = 25
A = b*h
Iy = b*h^3/12
Iz = h*b^3/12
J = Iy+Iz
Kt = 6.5e-3
Ksy = Ksz = 5/6

# Expected stiffness matrix
S0 = diagm([E*A,G*A*Ksy,G*A*Ksz,G*J*Kt,E*Iy,E*Iz]);