using BeamROM

# Computes the generalized displacements of a cantilevered, prismatic, anisotropic beam subjected to end loads, using Hodges' Eq. (5.62)

# Beam properties
L = 0.55
EA = 9.59482e6
GA2 = 3e6
GA3 = 3e6
GJ = 6.13689
EI2 = 4.06643
EI3 = 2874.31

# Load case
loadCase = "all"

# Sectional stiffness matrix
S = anisotropic_stiffness_matrix(EA=EA,GA2=GA2,GA3=GA3,GJ=GJ,EI2=EI2,EI3=EI3)

# Generalized displacements
uθ = displacements_from_loads(loadCase=loadCase,L=L,S=S)

println("Generalized displacements:")
display(uθ)