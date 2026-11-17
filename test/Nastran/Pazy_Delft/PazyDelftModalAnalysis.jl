using BeamROM, LinearAlgebra, Printf
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

# Delft Pazy wing (nastran/Pazy_Delft/pazy_delft.bdf) modal (SOL 103) analysis with NASTRAN, with the tip rod rigid and flexible (the rod, property 13, tied to the spar tip). The .bdf of both are generated, and, if the results of one are not there yet, the other is still processed, so that all the .bdf to run are listed together; run again after copying the .op2 back (PazyDelft modal results: pazy_delft_modal_rigidtip.op2 and pazy_delft_modal_flextip.op2).

nd = 15
nModes = 15

tips = (:rigid, :flexible)
responses = pazy_tip_responses(tip -> nastran_pazy_delft_modal_response(; nd=nd, tip=tip), tips)
λ, _ = responses[:rigid]
λFlex, _ = responses[:flexible]

for tip in tips
    foreach(file -> println("Saved ", file), nastran_pazy_delft_modal_shapes(; tip=tip, nmodes=nModes))
end

nd >= nModes || error("nd = $nd must be at least $nModes")
f = sqrt.(λ[1:nModes]) / (2π)
fFlex = sqrt.(λFlex[1:nModes]) / (2π)
println("\nNatural frequencies of the Delft Pazy wing (Nastran) with the rigid and with the flexible tip rod:")
@printf("%-4s %14s %14s %10s\n", "Mode", "Rigid [Hz]", "Flexible [Hz]", "Diff [%]")
for i in 1:nModes
    @printf("%-4d %14.4f %14.4f %10.2f\n", i, f[i], fFlex[i], 100 * (fFlex[i] - f[i]) / f[i])
end
