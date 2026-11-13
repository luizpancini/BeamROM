using BeamROM, LinearAlgebra, Printf, PrettyTables
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

# Equivalent beam inertia model of the Technion Pazy wing Nastran FEM, condensed to a reference axis.

# Sweep angle [deg] of the wing (S0 is the straight wing)
sweep = 0

# Position of the reference axis, as a fraction of the chord
reference_axis = 0.4653

# Beam mesh: :nominal places the nodes at the root, at the ribs (the first 2 ribs, close to the root, are lumped in the root node, and the last 3 of the 4 close to the tip are lumped in the tip node) and at the tip of the spar. :refined adds a node at the middle of each segment of :nominal, except the last one (twice as many elements, minus one). An integer gives that number of equally spaced nodes (e.g. 16) and a vector gives fractions of the span
# Each combination of reference axis and mesh is saved in its own file (S0_RA<reference_axis>_inertia_matrices_<nominal|refined>mesh.csv)
meshes = [:nominal, :refined]

for mesh in meshes

    # Get GFEM model and compute the equivalent beam mass model
    bdf = joinpath(pkgdir(BeamROM), "nastran", "Pazy_Technion", "S0", "S0.bdf")
    model = pazy_beam_mass_model(bdf; nodes=mesh, reference_axis, sweep, verbose=false)

    # Beam node coordinates in the output: x and z are in the Nastran basic frame (x = chordwise position from the origin of the reference system, i.e. the leading edge at x = 0 plus reference_axis of the chord), and only y starts at zero at the first node (the root of the spar axis). The offsets and tensors are relative to the nodes and to the centers of mass
    xref0 = [[x[1], x[2] - model.xref[1][2], x[3]] for x in model.xref]
    write_beam_mass_csv(pazy_mass_model_file(; reference_axis, mesh), xref0, model.m, model.offsets, model.tensors)

    # Print summary of the reference axis and mass properties
    total_mass = sum(model.m)
    nodes = length(model.xref)
    @printf("Reference axis in the Nastran frame: x = %.5f m, z = %.5f m, y from %.5f to %.5f m (%d nodes); output y relative to the first node\n", model.xref[1][1], model.xref[1][3], model.xref[1][2], model.xref[end][2], nodes)
    @printf("Total mass: beam model %.8f kg, lumped GFEM %.8f kg (difference %.2e %%)\n", total_mass, model.total_mass, 100 * (total_mass - model.total_mass) / model.total_mass)

    # Per node: y position [m], mass [g], center of mass offset from the reference axis [mm] and inertia tensor components about the center of mass [kg mm²]. Jyy is the torsional inertia about the span axis
    table = [1:nodes [x[2] for x in xref0] 1e3 * model.m 1e3 * [o[1] for o in model.offsets] 1e3 * [o[2] for o in model.offsets] 1e3 * [o[3] for o in model.offsets] 1e6 * [t[1, 1] for t in model.tensors] 1e6 * [t[2, 2] for t in model.tensors] 1e6 * [t[3, 3] for t in model.tensors]]
    pretty_table(table;
                title="Technion Pazy wing equivalent beam mass model",
                fit_table_in_display_horizontally=false,
                fit_table_in_display_vertically=false,
                column_labels=["Node", "y [m]", "m [g]", "Δx [mm]", "Δy [mm]", "Δz [mm]", "Jxx [kg mm²]", "Jyy [kg mm²]", "Jzz [kg mm²]"],
                formatters=[fmt__printf("%d", [1]), fmt__printf("%.5f", [2]), fmt__printf("%.3f", [3]), fmt__printf("%.3f", [4, 5, 6]), fmt__printf("%.4e", [7, 8, 9])])
end