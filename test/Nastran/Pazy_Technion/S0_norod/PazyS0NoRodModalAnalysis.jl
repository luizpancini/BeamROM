using BeamROM, LinearAlgebra, Printf
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

# Technion Pazy wing without the tip rod, modal (SOL 103) analysis with NASTRAN. The natural frequencies of the first modes of the wing clamped at the root are printed and compared with those of a cantilever beam with the properties of the spar and no tip mass (the rod and its concentrated masses are not in the model). The frequencies of the modes are saved as S0_norod_modal_frequencies.csv (mode number, eigenvalue [(rad/s)²], frequency [Hz]) and the mode shapes as S0_norod_modal_shapes_wing.csv (translations of the sections of the spar, see extractModeShapes.py), in nastran/Pazy_Technion/S0_norod, to be copied to AeroBeams (test/referenceData/Pazy)

# Number of modes requested in the modal analysis (at least the number of modes to compare)
nd = 10
nModes = 10
nd >= nModes || error("nd = $nd must be at least $nModes")

# Eigenvalues [(rad/s)²] of the Nastran modal analysis and beam length
λ, L = nastran_pazy_no_rod_modal_response(; nd=nd)

# Zeros of f in (0, xmax): sign changes on a grid, refined by bisection
function zeros_in(f, xmax; dx=1e-3)
    xs = dx:dx:xmax
    fs = f.(xs)
    roots = Float64[]
    for i in 1:length(xs)-1
        if fs[i] * fs[i+1] < 0
            a, b = xs[i], xs[i+1]
            for _ in 1:60
                m = (a + b) / 2
                f(a) * f(m) <= 0 ? (b = m) : (a = m)
            end
            push!(roots, (a + b) / 2)
        end
    end
    return roots
end

# Determinant of the free-end conditions of a clamped-free beam without tip mass and tip rotary inertia, at λ = βL: its zeros satisfy cos(λ) cosh(λ) = -1
bending_det(λ) = 1 + cos(λ) * cosh(λ)

# Modes of the clamped-free spar without tip body: bending out of plane (Iy) and in plane (Iz), and torsion with x = ωL sqrt(ρJ / (GJKt)), whose free end has no inertia so that cos(x) = 0
function cantilever_modes(; L, E, G, ρ, A, Iy, Iz, J, Kt, nModes)
    xmax = (nModes + 1) * π
    λBending = zeros_in(bending_det, xmax)
    xTor = zeros_in(cos, xmax)
    min(length(λBending), length(xTor)) >= nModes || error("not enough modes below $xmax")
    modes = Tuple{String,Float64}[]
    for n in 1:nModes
        push!(modes, ("bending out of plane $n", λBending[n]^2 * sqrt(E * Iy / (ρ * A)) / L^2 / (2π)))
        push!(modes, ("bending in plane $n", λBending[n]^2 * sqrt(E * Iz / (ρ * A)) / L^2 / (2π)))
        push!(modes, ("torsion $n", xTor[n] / L * sqrt(G * J * Kt / (ρ * J)) / (2π)))
    end
    return first(sort(modes, by=last), nModes)
end

# Frequencies [Hz] of the expected cantilevered spar modes
b = 60e-3
h = 2.25e-3
E = 71e9
ν = 0.33
G = E / (2 * (1 + ν))
ρ = 2795
A = b*h
Iy = b*h^3/12
Iz = h*b^3/12
J = Iy+Iz
Kt = 6.5e-3
expectedModes = cantilever_modes(; L, E, G, ρ, A, Iy, Iz, J, Kt, nModes)

# Frequencies [Hz] of the first modes
f = sqrt.(λ[1:nModes]) / (2π)
println("\nNatural frequencies of the Pazy wing S0 without the tip rod (Nastran) and of the equivalent cantilever spar:")
@printf("%-4s %14s   %-26s %14s %10s\n", "Mode", "Nastran [Hz]", "Spar mode", "Spar [Hz]", "Diff [%]")
for (i, (name, fExpected)) in enumerate(expectedModes)
    @printf("%-4d %14.4f   %-26s %14.4f %10.2f\n", i, f[i], name, fExpected, 100 * (f[i] - fExpected) / fExpected)
end

# Frequencies and mode shapes of the first modes, to be copied to AeroBeams (test/referenceData/Pazy)
foreach(file -> println("Saved ", file), nastran_pazy_no_rod_modal_files(; nmodes=nModes))
