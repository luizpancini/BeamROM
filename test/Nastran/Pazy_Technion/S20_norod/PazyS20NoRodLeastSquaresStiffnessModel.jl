using BeamROM, LinearAlgebra, Printf, PrettyTables
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

# Technion Pazy wing with 20 degrees of sweep (S20) and no tip rod stiffness matrices by least squares (see PazyS0LeastSquaresStiffnessModel.jl): the compliance matrix of each segment between consecutive beam nodes (or of the whole beam) is fitted to the generalized displacements of all the beam nodes under the six unit loads at the end of the reference axis (the static analysis of PazyS20NoRodGlobalStiffnessModel.jl), subject to being positive semidefinite and bounded below (see least_squares_stiffness_matrices). The displacements of a node are the rigid body motion of the section of constant y (of the Nastran model) of the spar, referred to the point of the inclined reference axis there, in the beam frame rotated by the sweep angle, and the node is at the distance (y - y_root)/cos(sweep) from the root along the axis

# Sweep angle [deg]: the sections (planes of constant y) of the Nastran model are shifted along x by y tan(sweep), so the reference axis is inclined by this angle from the y axis
sweep = 20

# Version of the model (only the Baseline has the static analysis)
version = "Baseline"

# Position of the reference axis, as a fraction of the chord (the same along the span, so that the axis is inclined by the sweep angle).
reference_axis = 0.4565

# Beam mesh: :nominal places the nodes at the root, at the ribs (the first 2 ribs, close to the root, are lumped in the root node, and the last 3 of the 4 close to the tip are lumped in the tip node) and at the tip of the spar. :refined adds a node at the middle of each segment of :nominal, except the last one (twice as many elements, minus one). An integer gives that number of equally spaced nodes (e.g. 16) and a vector gives fractions of the span
# Each combination of reference axis, mesh and segments is saved in its own file (S20_<version>_norod_RA<reference_axis>_stiffness_matrices_<nominal|refined>mesh.csv, with the version in lower case and _global at the end for :global segments)
meshes = [:nominal, :refined]

# Segments: :nodes = one stiffness matrix for each segment between consecutive beam nodes; :global = a single one for the whole beam, fitted to all the nodes
segments = [:nodes, :global]

# Spar (60 x 2.25 mm aluminum plate): its shear stiffness G A (without correction factor) bounds the transverse shear stiffnesses, which the data do not determine (without a bound the compliance is singular, i.e. rigid shear)
b = 60e-3
h = 2.25e-3
E = 71e9
ν = 0.33
G = E / (2 * (1 + ν))
GAmax = G * b * h
Cmin = Diagonal([0, 1/GAmax, 1/GAmax, 0, 0, 0])

# Weight of the smoothness of the compliances along the beam (0 = none)
λ = 0

for mesh in meshes
    for segments in segments
        # Beam nodes (y of the Nastran model and coordinate along the inclined axis) and their generalized displacements under the six load cases
        bdf = pazy_swept_no_rod_bdf(sweep, version)
        xref = spar_reference_axis(read_gfem(bdf), mesh; reference_axis, sweep)
        y = [x[2] for x in xref]
        x1 = (y .- y[1]) ./ cosd(sweep)
        T, L = nastran_pazy_swept_no_rod_station_response(sweep, version, y; reference_axis)

        # Least squares stiffness matrices
        breaks = segments == :nodes ? x1 : [x1[1], x1[end]]
        fit = least_squares_stiffness_matrices(y=x1, T=T, breaks=breaks, Cmin=Cmin, λ=λ, details=true)
        S = fit.S

        # Save the stiffness matrices next to the Nastran files, with the coordinates of the ends of the segments along the axis, from the root
        write_discretized_stiffness_csv(pazy_swept_no_rod_least_squares_stiffness_file(sweep, version; reference_axis, mesh, segments), breaks, S)

        # Print summary
        @printf("S20 %s least squares stiffness matrices (%s): length of the beam axis L = %.6f m, %d nodes, %d segments, %d iterations\n", version, segments, L, length(y), length(S), fit.iterations)
        @printf("RMS of the weighted residuals of the generalized displacements: %.3e\n", fit.residual)
        @printf("All stiffness matrices positive definite: %s\n", all(isposdef(Symmetric(Sk)) for Sk in S))
        println("Generalized displacements of the last node: fitted / Nastran")
        display(round.(fit.T_fit[end] ./ T[end], digits=4))

        # Diagonal of the stiffness matrix of each segment, and the inverse of the diagonal of its compliance (the stiffness without the couplings)
        table = [1:length(S) breaks[1:end-1] breaks[2:end] reduce(hcat, diag(Sk) for Sk in S)' reduce(hcat, 1 ./ diag(Ck) for Ck in fit.C)']
        pretty_table(table;
                    title="Technion Pazy wing S20 ($version) least squares stiffness: diagonal of the stiffness matrix of each segment and inverse of the diagonal of its compliance",
                    fit_table_in_display_horizontally=false,
                    fit_table_in_display_vertically=false,
                    column_labels=["Segment", "Start x1", "End x1", "S11", "S22", "S33", "S44", "S55", "S66", "1/C11", "1/C22", "1/C33", "1/C44", "1/C55", "1/C66"],
                    formatters=[fmt__printf("%d", [1]), fmt__printf("%.5f", [2, 3]), fmt__printf("%.2e", 4:15)])
    end
end
