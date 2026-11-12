using BeamROM, LinearAlgebra, Printf, PrettyTables

# Equivalent beam mass model of the Nastran FEMs of the prismatic beam (CHEXA, CQUAD4 and CBAR models in /nastran), condensed to a user-specified reference axis (the elastic axis, which coincides with the shear center of this homogeneous rectangular section) and compared with the analytic mass properties of each spanwise bin

# Beam properties
L = 0.55
b = 60e-3
h = 2.25e-3
ρ = 2.795e3

# Number of beam nodes
nodes = 11

# Cases (element type, mesh ID and number of elements along the length, see nastran/generateBeam.py). The grids are assigned wholly to the nearest beam node, so the bin boundaries are only resolved up to one grid spacing L/nL
cases = [("solid", 1, 100), 
         ("solid", 2, 100), 
         ("solid", 3, 100), 
         ("solid", 4, 200), 
         ("shell", 1, 50), 
         ("line", 1, 10), 
         ("line", 4, 100)]

# Reference axis (y, z) in the Nastran basic frame: the beam is along x, with width along y and thickness along z. The solid has its section at y ∈ [0,b], z ∈ [0,h]; the shell mid-surface is at z = 0 (y ∈ [0,b]); the CBAR axis passes through the origin of the section
reference_axis = Dict("solid" => (b / 2, h / 2), 
                      "shell" => (b / 2, 0.0), 
                      "line" => (0.0, 0.0))

# Analytic mass properties of the bins centered at the beam nodes (half bins at the ends)
ΔL = L / (nodes - 1)
binwidth = [k in (1, nodes) ? 
            ΔL / 2 : 
            ΔL for k in 1:nodes]
bincenter = [k == 1 ? 
             ΔL / 4 : 
             k == nodes ? 
             L - ΔL / 4 : 
             (k - 1) * ΔL for k in 1:nodes]
m_exact = ρ * b * h .* binwidth
# Torsional inertia per unit mass of the rectangular section about its centroid
Jxx_per_mass = (b^2 + h^2) / 12

# Loop over each case and compare the condensed beam mass model with the analytic mass properties
for (element, mesh, nL) in cases
    yea, zea = reference_axis[element]
    xref = straight_axis([0.0, yea, zea], [L, yea, zea], nodes)
    bdf = joinpath(pkgdir(BeamROM), "nastran", "simple_beams", "$(element)_beam", "mesh$(mesh)", "$(element)_beam_static.bdf")
    model = beam_mass_model(bdf, xref; verbose=false)
    write_beam_mass_csv(joinpath(dirname(bdf), "$(element)_beam_mass_model.csv"), model.xref, model.m, model.offsets, model.tensors)

    err_Jxx = [abs(t[1, 1] / m - Jxx_per_mass) / Jxx_per_mass * 100 for (t, m) in zip(model.tensors, model.m)]
    # The section centroid lies on the reference axis
    err_yz = [norm(o[2:3]) / b * 100 for o in model.offsets]

    total_mass, total_exact = sum(model.m), ρ * L * b * h
    title = @sprintf("%s beam, mesh %d: total mass %.8f kg (exact %.8f kg, difference %.2e %%)", element, mesh, total_mass, total_exact, 100 * (total_mass - total_exact) / total_exact)

    # Per node: mass [g] and its percentage difference from the analytic bin mass, x offset of the center of mass [mm] (numerical and analytic), y-z offset as a percentage of the width b, and Jxx/m error [%]. The line model carries no section inertia, and the shell none through the thickness: only the solid is expected to match Jxx/m closely
    table = [1:nodes [x[1] for x in xref] 1e3 * model.m 1e3 * m_exact 100 * (model.m - m_exact) ./ m_exact 1e3 * [o[1] for o in model.offsets] 1e3 * (bincenter - [x[1] for x in xref]) err_yz err_Jxx]
    pretty_table(table;
                 title=title,
                 fit_table_in_display_horizontally=false,
                 fit_table_in_display_vertically=false,
                 column_labels=["Node", "x [m]", "m [g]", "m exact [g]", "Δm [%]", "Δx [mm]", "Δx exact [mm]", "Δyz/b [%]", "Jxx/m error [%]"],
                 formatters=[fmt__printf("%d", [1]), fmt__printf("%.4f", [2]), fmt__printf("%.4f", [3, 4]), fmt__printf("%.3f", [5]), fmt__printf("%.4f", [6, 7]), fmt__printf("%.1e", [8]), fmt__printf("%.3f", [9])])
end
