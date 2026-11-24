using BeamROM, LinearAlgebra, Printf, PrettyTables
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

# Pazy wing (Technion Nastran FEM, see nastran/Pazy_Technion/S0) stiffness matrices by least squares: the compliance matrix of each segment between consecutive beam nodes (or of the whole beam) is fitted to the generalized displacements of all the beam nodes under the six unit loads at the end of the reference axis (the static analysis of PazyS0GlobalStiffnessModel.jl), subject to being positive semidefinite and bounded below (see least_squares_stiffness_matrices). Unlike the exact inversion of PazyS0GlobalStiffnessModel.jl and PazyS0DiscretizedStiffnessModel.jl, this gives positive definite stiffness matrices: the transverse shear compliance of the slender wing is a tiny difference of large bending terms, dominated by the deviations from a prismatic beam (clamp, tip), which gives negative shear stiffnesses

# Sweep angle [deg] of the wing (S0 is the straight wing)
sweep = 0

# Position of the reference axis, as a fraction of the chord
reference_axis = 0.4653

# Beam mesh: :nominal places the nodes at the root, at the ribs (the first 2 ribs, close to the root, are lumped in the root node, and the last 3 of the 4 close to the tip are lumped in the tip node) and at the tip of the spar. :refined adds a node at the middle of each segment of :nominal, except the last one (twice as many elements, minus one). An integer gives that number of equally spaced nodes (e.g. 16) and a vector gives fractions of the span
# Each combination of reference axis, mesh, tip and segments is saved in its own file (S0_RA<reference_axis>_stiffness_matrices_<nominal|refined>mesh_<flex|rigid>tip.csv, with _global at the end for :global segments)
meshes = [:nominal, :refined]

# Tip rod of the Nastran model: :rigid (the rod, property 13, is made of a rigid nylon, so the tip section of the wing is, as the stiffness matrices require) or :flexible (the rod of the model, which bends). The results of each have the suffix _rigidtip or _flextip in their file names (see pazy_tip_suffix)
tips = [:flexible, :rigid]

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
    for tip in tips
        for segments in segments
            # Beam nodes and their generalized displacements under the six load cases
            bdf = joinpath(pkgdir(BeamROM), "nastran", "Pazy_Technion", "S0", "S0.bdf")
            xref = spar_reference_axis(read_gfem(bdf), mesh; reference_axis, sweep)
            y = [x[2] for x in xref]
            y0 = y .- y[1]
            T, L = nastran_pazy_station_response(y; reference_axis, tip)

            # Least squares stiffness matrices
            breaks = segments == :nodes ? y0 : [y0[1], y0[end]]
            fit = least_squares_stiffness_matrices(y=y0, T=T, breaks=breaks, Cmin=Cmin, λ=λ, details=true)
            S = fit.S

            # Save the stiffness matrices next to the Nastran files, with the coordinates of the ends of the segments relative to the first node (the root of the spar axis)
            write_discretized_stiffness_csv(pazy_least_squares_stiffness_file(; reference_axis, mesh, segments, tip), breaks, S)

            # Print summary
            @printf("Least squares stiffness matrices (%s): %d nodes, %d segments, %d iterations\n", segments, length(y), length(S), fit.iterations)
            @printf("RMS of the weighted residuals of the generalized displacements: %.3e\n", fit.residual)
            @printf("All stiffness matrices positive definite: %s\n", all(isposdef(Symmetric(Sk)) for Sk in S))
            println("Generalized displacements of the last node: fitted / Nastran")
            display(round.(fit.T_fit[end] ./ T[end], digits=4))

            # Diagonal of the stiffness matrix of each segment, and the inverse of the diagonal of its compliance (the stiffness without the couplings)
            table = [1:length(S) breaks[1:end-1] breaks[2:end] reduce(hcat, diag(Sk) for Sk in S)' reduce(hcat, 1 ./ diag(Ck) for Ck in fit.C)']
            pretty_table(table;
                        title="Technion Pazy wing least squares stiffness: diagonal of the stiffness matrix of each segment and inverse of the diagonal of its compliance",
                        fit_table_in_display_horizontally=false,
                        fit_table_in_display_vertically=false,
                        column_labels=["Segment", "y1", "y2", "S11", "S22", "S33", "S44", "S55", "S66", "1/C11", "1/C22", "1/C33", "1/C44", "1/C55", "1/C66"],
                        formatters=[fmt__printf("%d", [1]), fmt__printf("%.5f", [2, 3]), fmt__printf("%.2e", 4:15)])
        end
    end
end
