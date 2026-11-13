using BeamROM, LinearAlgebra, Printf, PrettyTables
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

sweep = 0
reference_axis = 0.445

# Beam mesh: :nominal places the nodes at the root, at the first 14 ribs (the last 3 of the 4 ribs close to the tip are lumped in the tip node) and at the tip of the spar, for 15 elements, the last one from the 14th rib to the tip. :refined adds a node at the middle of each element of :nominal, except the last one (twice as many elements, minus one)
meshes = [:nominal, :refined]

for mesh in meshes

    model = pazy_delft_beam_mass_model(; nodes=mesh, reference_axis, sweep, verbose=false)

    xref0 = [[x[1], x[2] - model.xref[1][2], x[3]] for x in model.xref]
    write_beam_mass_csv(pazy_delft_mass_model_file(; reference_axis, mesh), xref0, model.m, model.offsets, model.tensors)

    # Print summary of the reference axis and mass properties
    total_mass = sum(model.m)
    nodes = length(model.xref)
    @printf("Reference axis in the Nastran frame: x = %.5f m, z = %.5f m, y from %.5f to %.5f m (%d nodes)\n", model.xref[1][1], model.xref[1][3], model.xref[1][2], model.xref[end][2], nodes)
    @printf("Total mass: beam model %.8f kg, lumped GFEM %.8f kg (difference %.2e %%)\n", total_mass, model.total_mass, 100 * (total_mass - model.total_mass) / model.total_mass)

    table = [1:nodes [x[2] for x in xref0] 1e3 * model.m 1e3 * [o[1] for o in model.offsets] 1e3 * [o[2] for o in model.offsets] 1e3 * [o[3] for o in model.offsets] 1e6 * [t[1, 1] for t in model.tensors] 1e6 * [t[2, 2] for t in model.tensors] 1e6 * [t[3, 3] for t in model.tensors]]
    pretty_table(table;
                title="Delft Pazy wing equivalent beam mass model",
                fit_table_in_display_horizontally=false,
                fit_table_in_display_vertically=false,
                column_labels=["Node", "y [m]", "m [g]", "Δx [mm]", "Δy [mm]", "Δz [mm]", "Jxx [kg mm²]", "Jyy [kg mm²]", "Jzz [kg mm²]"],
                formatters=[fmt__printf("%d", [1]), fmt__printf("%.5f", [2]), fmt__printf("%.3f", [3]), fmt__printf("%.3f", [4, 5, 6]), fmt__printf("%.4e", [7, 8, 9])])
end
