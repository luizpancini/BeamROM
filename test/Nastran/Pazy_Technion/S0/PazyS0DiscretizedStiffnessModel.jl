using BeamROM, LinearAlgebra, Printf, PrettyTables
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

# Pazy wing (Technion Nastran FEM see nastran/Pazy_Technion/S0) discretized stiffness condensation: instead of one stiffness matrix for the whole beam (see PazyS0StiffnessModel.jl), one for each segment between consecutive beam nodes. The static analysis (SOL 101) is the same: the wing is clamped at the root section and the six unit loads are applied at the end of the reference axis. The generalized displacements of every beam node (not only of the last one) are taken under the six load cases, and each segment is treated as a cantilever clamped at its upstream node, loaded at its downstream node by the generalized forces transmitted through it, with the displacements of the downstream node relative to the upstream one (in the frame of the rigid body rotation of the upstream node)

# Sweep angle [deg] of the wing (S0 is the straight wing)
sweep = 0

# Position of the reference axis, as a fraction of the chord
reference_axis = 0.4653

# Beam mesh: :nominal places the nodes at the root, at the ribs (the first 2 ribs, close to the root, are lumped in the root node, and the last 3 of the 4 close to the tip are lumped in the tip node) and at the tip of the spar. :refined adds a node at the middle of each segment of :nominal, except the last one (twice as many elements, minus one). An integer gives that number of equally spaced nodes (e.g. 16) and a vector gives fractions of the span
# Each combination of reference axis, mesh and tip is saved in its own file (S0_RA<reference_axis>_discretized_stiffness_matrices_<nominal|refined>mesh_<flex|rigid>tip.csv)
mesh = :refined

# Tip rod of the Nastran model: :rigid (the rod is rigid, so the tip section of the wing is, as the stiffness matrices require) or :flexible (the rod of the model, which bends). The results of each have the suffix _rigidtip or _flextip in their file names (see pazy_tip_suffix)
tip = :rigid

# Get the beam nodes of the mass model and the generalized displacements of each one under the six load cases
bdf = joinpath(pkgdir(BeamROM), "nastran", "Pazy_Technion", "S0", "S0.bdf")
xref = spar_reference_axis(read_gfem(bdf), mesh; reference_axis, sweep)
y = [x[2] for x in xref]
nodes = length(y)
T, L = nastran_pazy_station_response(y; reference_axis, tip)

# Compute the stiffness matrix of each segment
S = discretized_stiffness_matrices(y=y, T=T)

# Save the stiffness matrices next to the Nastran files, with the beam node coordinates relative to the first node (the root of the spar axis). They are not symmetric, so all their entries are saved
y0 = y .- y[1]
write_discretized_stiffness_csv(pazy_discretized_stiffness_file(; reference_axis, mesh, tip), y0, S; symmetric=false)

# Stiffness matrix of the whole beam, from the generalized displacements of the last node, and of the spar alone (60 x 2.25 mm aluminum plate, see OneDBeamNastran.jl for the correction factors), as references
Tend = T[end]
Sglobal = stiffness_matrix(L=L,u_F1=Tend[1:3,1],u_F2=Tend[1:3,2],u_F3=Tend[1:3,3],u_M1=Tend[1:3,4],u_M2=Tend[1:3,5],u_M3=Tend[1:3,6],θ_F1=Tend[4:6,1],θ_F2=Tend[4:6,2],θ_F3=Tend[4:6,3],θ_M1=Tend[4:6,4],θ_M2=Tend[4:6,5],θ_M3=Tend[4:6,6])
b = 60e-3
h = 2.25e-3
E = 71e9
ν = 0.33
G = E / (2 * (1 + ν))
A = b*h
Iy = b*h^3/12
Iz = h*b^3/12
J = Iy+Iz
Kt = 6.5e-3
Ksy = Ksz = 5/6
S0 = diagm([E*A,G*A*Ksy,G*A*Ksz,G*J*Kt,E*Iy,E*Iz])

# Print summary
@printf("Generalized displacements of the sections (%d nodes, %d segments)\n", nodes, nodes-1)
# Asymmetry of the segment stiffness matrices: zero for a prismatic beam, nonzero when the nodal displacements are not exactly reciprocal (e.g. with warping)
@printf("Largest asymmetry of the segment stiffness matrices: %.2e (relative)\n", maximum(norm(Sk - Sk') / norm(Sk) for Sk in S))

# Diagonal of the stiffness matrix of each segment (EA, GA2, GA3, GJ, EI2, EI3, as in the spar stiffness matrix), y positions relative to the first node [m] and length [mm]
diags = reduce(hcat, diag(Sk) for Sk in S)'
table = [1:nodes-1 y0[1:end-1] y0[2:end] 1e3 * diff(y) diags]
pretty_table(table;
             title="Technion Pazy wing discretized stiffness: diagonal of each segment",
             fit_table_in_display_horizontally=false,
             fit_table_in_display_vertically=false,
             column_labels=["Segment", "y1 [m]", "y2 [m]", "Length [mm]", "S11 [N]", "S22 [N]", "S33 [N]", "S44 [N m²]", "S55 [N m²]", "S66 [N m²]"],
             formatters=[fmt__printf("%d", [1]), fmt__printf("%.5f", [2, 3]), fmt__printf("%.2f", [4]), fmt__printf("%.4e", 5:10)])

# Ratio of the diagonal of each segment to that of the spar alone, and of the whole beam
ratios = [diags ./ diag(S0)'; (diag(Sglobal) ./ diag(S0))']
pretty_table(ratios;
             title="Diagonal of the stiffness matrices relative to the spar alone (last row: whole beam)",
             fit_table_in_display_horizontally=false,
             fit_table_in_display_vertically=false,
             row_labels=[string.(1:nodes-1); "Beam"],
             column_labels=["S11", "S22", "S33", "S44", "S55", "S66"],
             formatters=[fmt__printf("%.3f", 1:6)])
