using BeamROM

# Computes the generalized displacements of a cantilevered, prismatic, anisotropic beam subjected to end loads, using Hodges' Eq. (5.62)

# Beam properties
L = 1
EA = 2.00136e7
GA2 = 7.35111e6
GA3 = 4.17224e4
GJ = 10.0993
EI2 = 6.71671
EI3 = 4168.17

# Load case
loadCase = "all"

# Sectional stiffness matrix
S = anisotropic_stiffness_matrix(EA=EA,GA2=GA2,GA3=GA3,GJ=GJ,EI2=EI2,EI3=EI3)

# Generalized displacements
uθ = displacements_from_loads(loadCase=loadCase,L=L,S=S)

println("Generalized displacements:")
display(uθ)