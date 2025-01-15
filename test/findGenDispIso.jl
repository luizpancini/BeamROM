using BeamROM

# Computes the generalized displacements of a cantilevered, prismatic, isotropic beam subjected to end loads, using Hodges' Eq. (5.62)

# Beam properties
L = 1
EA = 2e7
GA2 = 6.4e6
GA3 = GA2
GJ = 10.2
EI2 = 6.666667
EI3 = 4166.666667

# Load case
loadCase = "all"

# Sectional stiffness matrix
S = isotropic_stiffness_matrix(EA=EA,GA2=GA2,GA3=GA3,GJ=GJ,EI2=EI2,EI3=EI3)

# Generalized displacements
uθ = displacements_from_loads(loadCase=loadCase,L=L,S=S)

println("Generalized displacements:")
display(uθ)