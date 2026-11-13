using BeamROM, LinearAlgebra, Printf, PrettyTables
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

# Equivalent beam inertia model of the Technion Pazy wing with 20 degrees of sweep (S20) and no tip rod Nastran FEM, condensed to a reference axis. The model has three versions, with the tip weights at the leading edge (LE) or at the trailing edge (TE), and the Baseline, all of them without the tip rod

# Versions of the model ("LE", "TE" and/or "Baseline")
versions = ["LE", "TE", "Baseline"]

# Sweep angle [deg]: the sections (planes of constant y) of the Nastran model are shifted along x by y tan(sweep), so the reference axis is inclined by this angle from the y axis
sweep = 20

# Position of the reference axis, as a fraction of the chord (the same along the span, so that the axis is inclined by the sweep angle).
reference_axis = 0.4565

# Beam node distribution: :ribs places the nodes at the root, at the ribs (the first 2 ribs, close to the root, are lumped in the root node, and the last 3 of the 4 close to the tip are lumped in the tip node) and at the tip of the spar. :ribs_halved adds a node at the middle of each segment of :ribs, except the last one (twice as many elements, minus one). An integer gives that number of equally spaced nodes (e.g. 16) and a vector gives fractions of the span
# Each combination of reference_axis and nodes_distribution is saved in its own file, in the NoRod folder (S20_<version>_norod_RA<reference_axis>_ND<nodes_distribution>_beam_mass_model.csv)
nodes_distribution = [:ribs, :ribs_halved]

for version in versions
    for mesh in nodes_distribution
        # Bulk data of the model without the tip rod (see pazy_swept_no_rod_bdf) and the equivalent beam mass model
        bdf = pazy_swept_no_rod_bdf(sweep, version)
        model = pazy_beam_mass_model(bdf; nodes=mesh, reference_axis, sweep, verbose=false)

        # Beam node coordinates in the output: x and z are in the Nastran basic frame (x = chordwise position from the origin of the reference system, which varies along the span with the sweep), and only y starts at zero at the first node (the root of the spar axis). The offsets and tensors are relative to the nodes and to the centers of mass, in the Nastran frame (x chordwise, y spanwise, z normal), not in the frame of the inclined axis
        xref0 = [[x[1], x[2] - model.xref[1][2], x[3]] for x in model.xref]
        write_beam_mass_csv(pazy_swept_no_rod_mass_model_file(sweep, version; reference_axis, nodes_distribution=mesh), xref0, model.m, model.offsets, model.tensors)

        # Print summary of the reference axis and mass properties
        total_mass = sum(model.m)
        nodes = length(model.xref)
        @printf("S20 %s: reference axis in the Nastran frame, from x = %.5f m, y = %.5f m to x = %.5f m, y = %.5f m, at z = %.5f m (%d nodes, inclined by %.1f°); output y relative to the first node\n", version, model.xref[1][1], model.xref[1][2], model.xref[end][1], model.xref[end][2], model.xref[1][3], nodes, sweep)
        @printf("Total mass: beam model %.8f kg, lumped GFEM %.8f kg (difference %.2e %%)\n", total_mass, model.total_mass, 100 * (total_mass - model.total_mass) / model.total_mass)

        # Per node: x and y positions [m], mass [g], center of mass offset from the reference axis [mm] and inertia tensor components about the center of mass [kg mm²], in the Nastran frame. Jyy is the inertia about the y axis, which is not the span axis of the beam (the axis is inclined)
        table = [1:nodes [x[1] for x in xref0] [x[2] for x in xref0] 1e3 * model.m 1e3 * [o[1] for o in model.offsets] 1e3 * [o[2] for o in model.offsets] 1e3 * [o[3] for o in model.offsets] 1e6 * [t[1, 1] for t in model.tensors] 1e6 * [t[2, 2] for t in model.tensors] 1e6 * [t[3, 3] for t in model.tensors]]
        pretty_table(table;
                    title="No rod Technion Pazy wing S20 ($version) equivalent beam mass model",
                    fit_table_in_display_horizontally=false,
                    fit_table_in_display_vertically=false,
                    column_labels=["Node", "x [m]", "y [m]", "m [g]", "Δx [mm]", "Δy [mm]", "Δz [mm]", "Jxx [kg mm²]", "Jyy [kg mm²]", "Jzz [kg mm²]"],
                    formatters=[fmt__printf("%d", [1]), fmt__printf("%.5f", [2, 3]), fmt__printf("%.3f", [4]), fmt__printf("%.3f", [5, 6, 7]), fmt__printf("%.4e", [8, 9, 10])])
    end
end
