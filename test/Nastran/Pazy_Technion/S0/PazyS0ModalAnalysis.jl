using BeamROM, LinearAlgebra, Printf
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

# Technion Pazy wing modal (SOL 103) analysis with NASTRAN, with the tip rod rigid (the tip section of the wing is rigid, as the beam model of the stiffness matrices) and flexible (the rod of the model, which bends). The .bdf of both are generated, and, if the results of one are not there yet, the others are still processed, so that all the .bdf to run (e.g. on a Windows VM) are listed together; run again after copying the .op2 back (S0_modal_rigidtip.op2 and S0_modal_flextip.op2). The eigenvalues of the first modes of the wing with the rigid tip clamped at the root are compared with those of a cantilever beam with the properties of the spar and the tip mass (tip rod and concentrated masses), and those of the two tips with each other. The frequencies of each are saved as S0_modal_frequencies_<rigidtip or flextip>.csv (mode number, eigenvalue [(rad/s)²], frequency [Hz]), to be copied to AeroBeams (test/referenceData/Pazy), which compares those of the flexible tip with the beam model with the rod modeled as a beam (PazyWingBeamRodModal.jl), and the mode shapes as S0_modal_shapes_<rigidtip or flextip>_<wing or rod>.csv (translations of the sections of the spar and of the rod, see extractModeShapes.py), which it uses to pair the modes by their shapes (modal assurance criterion)

# Number of modes requested in the modal analysis (at least the number of modes to compare)
nd = 15
nModes = 15

# Eigenvalues [(rad/s)²] of the Nastran modal analysis with the rigid and the flexible tip rod, and beam length
tips = (:rigid, :flexible)
responses = pazy_tip_responses(tip -> nastran_pazy_modal_response(; nd=nd, tip=tip), tips)
λ, L = responses[:rigid] # the one of the cantilever spar with the tip body below
λFlex, _ = responses[:flexible]

# Mode shapes of the first modes with each tip, saved as S0_modal_shapes_<rigidtip or flextip>_<wing or rod>.csv, to be copied to AeroBeams (test/referenceData/Pazy)
for tip in tips
    foreach(file -> println("Saved ", file), nastran_pazy_modal_shapes(; tip=tip, nmodes=nModes))
end

# # Zeros of f in (0, xmax): sign changes on a grid, refined by bisection
# function zeros_in(f, xmax; dx=1e-3)
#     xs = dx:dx:xmax
#     fs = f.(xs)
#     roots = Float64[]
#     for i in 1:length(xs)-1
#         if fs[i] * fs[i+1] < 0
#             a, b = xs[i], xs[i+1]
#             for _ in 1:500
#                 m = (a + b) / 2
#                 f(a) * f(m) <= 0 ? (b = m) : (a = m)
#             end
#             push!(roots, (a + b) / 2)
#         end
#     end
#     return roots
# end

# # Determinant of the free-end conditions of a clamped-free beam with tip mass μ = m/(ρAL) and tip rotary inertia γ = I/(ρAL³), at λ = βL
# function bending_det(λ, μ, γ)
#     c, s, ch, sh = cos(λ), sin(λ), cosh(λ), sinh(λ)
#     a11 = -(c + ch) + γ * λ^3 * (s + sh)
#     a12 = -(s + sh) - γ * λ^3 * (c - ch)
#     a21 = (s - sh) + μ * λ * (c - ch)
#     a22 = -(c + ch) + μ * λ * (s - sh)
#     return a11 * a22 - a12 * a21
# end

# # Modes of the clamped-free spar with a tip body of mass `mTip` and rotary inertias `IxTip` (flap), `IyTip` (torsion) and `IzTip` (edge) about the elastic axis at the free end
# function cantilever_modes(; L, E, G, ρ, A, Iy, Iz, J, Kt, nModes, mTip, IxTip, IyTip, IzTip)
#     xmax = (nModes + 1) * π
#     μ = mTip / (ρ * A * L)
#     λOut = zeros_in(λ -> bending_det(λ, μ, IxTip / (ρ * A * L^3)), xmax)
#     λIn = zeros_in(λ -> bending_det(λ, μ, IzTip / (ρ * A * L^3)), xmax)
#     # Torsion with x = ωL sqrt(ρJ / (GJKt)): the tip inertia gives IyTip x sin(x) = ρJL cos(x)
#     xTor = zeros_in(x -> IyTip * x * sin(x) - ρ * J * L * cos(x), xmax)
#     min(length(λOut), length(λIn), length(xTor)) >= nModes || error("not enough modes below $xmax")
#     modes = Tuple{String,Float64}[]
#     for n in 1:nModes
#         push!(modes, ("bending out of plane $n", λOut[n]^2 * sqrt(E * Iy / (ρ * A)) / L^2 / (2π)))
#         push!(modes, ("bending in plane $n", λIn[n]^2 * sqrt(E * Iz / (ρ * A)) / L^2 / (2π)))
#         push!(modes, ("torsion $n", xTor[n] / L * sqrt(G * J * Kt / (ρ * J)) / (2π)))
#     end
#     return first(sort(modes, by=last), nModes)
# end

# # Frequencies [Hz] of the expected cantilevered spar modes
# b = 60e-3
# h = 2.25e-3
# E = 71e9
# ν = 0.33
# G = E / (2 * (1 + ν))
# ρ = 2795
# A = b*h
# Iy = b*h^3/12
# Iz = h*b^3/12
# J = Iy+Iz
# Kt = 6.5e-3
# # Tip body at the spar end (from S0.bdf): rod PBEAM 13 (solid circular section, tied rigidly to the spar end) and CONM2 masses 669 and 200669 at its center (node 1212)
# ρRod = 930
# RRod = 5e-3
# LRod = 0.2682
# mRod = ρRod * π * RRod^2 * LRod
# mTip = mRod + 0.015 + 0.0011
# # Offsets [m] of the rod center from the spar elastic axis at its free end, in chord and span
# dxTip = 0.05 - 0.043494
# dyTip = 0.5568333 - 0.5519937
# IxTip = mRod * RRod^2 / 2 + mTip * dyTip^2
# IyTip = mRod * (LRod^2 / 12 + RRod^2 / 4) + mTip * dxTip^2
# IzTip = mRod * (LRod^2 / 12 + RRod^2 / 4) + mTip * (dxTip^2 + dyTip^2)
# expectedModes = cantilever_modes(; L, E, G, ρ, A, Iy, Iz, J, Kt, nModes, mTip, IxTip, IyTip, IzTip)

# Natural frequencies of the wing with the rigid and with the flexible tip rod
nd >= nModes || error("nd = $nd must be at least $nModes")
f = sqrt.(λ[1:nModes]) / (2π)
fFlex = sqrt.(λFlex[1:nModes]) / (2π)
println("\nNatural frequencies of the Pazy wing (Nastran) with the rigid and with the flexible tip rod:")
@printf("%-4s %14s %14s %10s\n", "Mode", "Rigid [Hz]", "Flexible [Hz]", "Diff [%]")
for i in 1:nModes
    @printf("%-4d %14.4f %14.4f %10.2f\n", i, f[i], fFlex[i], 100 * (fFlex[i] - f[i]) / f[i])
end
