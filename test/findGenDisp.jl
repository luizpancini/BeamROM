using BeamROM

# Computes the generalized displacements of a cantilevered, prismatic, anisotropic beam subjected to end loads, using Hodges' Eq. (5.62)

# Beam properties
L = 1
EA = 2e7
GAy = 6.4e6
GAz = GAy
GJ = 10.2
EIy = 6.666667
EIz = 4166.666667

# Load case
loadCase = "all"

# Sectional stiffness matrix
S = isotropic_stiffness_matrix(EA=EA,GAy=GAy,GAz=GAz,GJ=GJ,EIy=EIy,EIz=EIz)
S[3,5] = S[5,3] = -0
S[1,4] = S[4,1] = 0

# Generalized displacements
uθ = displacements_from_loads(loadCase=loadCase,L=L,S=S)

println("Generalized displacements:")
display(uθ)