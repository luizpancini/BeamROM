using BeamROM, LinearAlgebra, Printf, PrettyTables
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

# Delft Pazy wing (nastran/Pazy_Delft/pazy_delft.bdf) stiffness matrices by least squares: the compliance matrix of each segment between consecutive beam nodes (or of the whole beam) is fitted to the generalized displacements of all the beam nodes under the six unit loads at the end of the reference axis (the static analysis of PazyDelftGlobalStiffnessModel.jl), subject to being positive semidefinite and bounded below (see least_squares_stiffness_matrices).

reference_axis = 0.445

# Beam mesh: :nominal places the nodes at the root, at the first 14 ribs (the last 3 of the 4 ribs close to the tip are lumped in the tip node) and at the tip of the spar, for 15 elements, the last one from the 14th rib to the tip. :refined adds a node at the middle of each element of :nominal, except the last one (twice as many elements, minus one)
meshes = [:nominal, :refined]

tips = [:flexible, :rigid]

segments = [:nodes, :global, :ends]

# Spar (60 x 1.8 mm aluminum plate): its shear stiffness G A (without correction factor) bounds the transverse shear stiffnesses, which the data do not determine (without a bound the compliance is singular, i.e. rigid shear)
b = 60e-3
h = 1.8e-3
E = 71e9
ν = 0.33
G = E / (2 * (1 + ν))
GAmax = G * b * h
Cmin = Diagonal([0, 1/GAmax, 1/GAmax, 0, 0, 0])

λ = 0

for mesh in meshes
    for tip in tips
        for segments in segments
            xref = pazy_delft_spar_reference_axis(mesh; reference_axis)
            y = [x[2] for x in xref]
            y0 = y .- y[1]
            T, L = nastran_pazy_delft_station_response(y; reference_axis, tip)

            breaks = if segments == :nodes
                y0
            elseif segments == :ends
                [y0[1], y0[2], y0[3], y0[end-2], y0[end-1], y0[end]]
            else
                [y0[1], y0[end]]
            end
            fit = least_squares_stiffness_matrices(y=y0, T=T, breaks=breaks, Cmin=Cmin, λ=λ, details=true)
            S = fit.S

            write_discretized_stiffness_csv(pazy_delft_least_squares_stiffness_file(; reference_axis, mesh, segments, tip), breaks, S)

            @printf("Least squares stiffness matrices (%s, %s tip): %d nodes, %d segments, %d iterations\n", segments, tip, length(y), length(S), fit.iterations)
            @printf("RMS of the weighted residuals of the generalized displacements: %.3e\n", fit.residual)
            @printf("All stiffness matrices positive definite: %s\n", all(isposdef(Symmetric(Sk)) for Sk in S))
            println("Generalized displacements of the last node: fitted / Nastran")
            display(round.(fit.T_fit[end] ./ T[end], digits=4))

            table = [1:length(S) breaks[1:end-1] breaks[2:end] reduce(hcat, diag(Sk) for Sk in S)' reduce(hcat, 1 ./ diag(Ck) for Ck in fit.C)']
            pretty_table(table;
                        title="Delft Pazy wing least squares stiffness ($tip tip): diagonal of the stiffness matrix of each segment and inverse of the diagonal of its compliance",
                        fit_table_in_display_horizontally=false,
                        fit_table_in_display_vertically=false,
                        column_labels=["Segment", "y1", "y2", "S11", "S22", "S33", "S44", "S55", "S66", "1/C11", "1/C22", "1/C33", "1/C44", "1/C55", "1/C66"],
                        formatters=[fmt__printf("%d", [1]), fmt__printf("%.5f", [2, 3]), fmt__printf("%.2e", 4:15)])
        end
    end
end
