using BeamROM, LinearAlgebra, Printf
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

# Technion Pazy wing with 10 degrees of sweep (S10, versions LE and TE of the tip weights, and Baseline) modal (SOL 103) analysis with NASTRAN, with the tip rod rigid (a rigid nylon in the properties 13, the rod, and 24, the breaking beam, so that the tip section of the wing is rigid, as the beam model of the stiffness matrices) and flexible (the materials of the model, so the rod and the breaking beam bend). The .bdf of every version and tip are generated, and, if the results of one are not there yet, the others are still processed, so that all the .bdf to run (e.g. on a Windows VM) are listed together; run again after copying the .op2 back (S10_<version>_modal_rigidtip.op2 and S10_<version>_modal_flextip.op2, in the folder of each version). The natural frequencies of the first modes of the wing clamped at the root with each tip are printed, and the frequencies of each are saved as S10_<version>_modal_frequencies_<rigidtip or flextip>.csv (mode number, eigenvalue [(rad/s)²], frequency [Hz]) and the mode shapes as S10_<version>_modal_shapes_<rigidtip or flextip>_<wing or rod>.csv (translations of the sections of the spar and of the rod, see extractModeShapes.py), to be copied to AeroBeams (test/referenceData/Pazy)

# Versions of the model ("LE", "TE" and/or "Baseline")
versions = ["LE", "TE", "Baseline"]

# Sweep angle [deg]: the sections (planes of constant y) of the Nastran model are shifted along x by y tan(sweep), so the axis of the spar is inclined by this angle from the y axis
sweep = 10

# Number of modes requested in the modal analysis (at least the number of modes to compare)
nd = 15
nModes = 15
nd >= nModes || error("nd = $nd must be at least $nModes")

# Eigenvalues [(rad/s)²] of the Nastran modal analysis, for each version and tip rod (rigid and flexible)
tips = (:rigid, :flexible)
responses = pazy_swept_tip_responses((version, tip) -> nastran_pazy_swept_modal_response(sweep, version; nd, tip), sweep, versions, tips)

for version in versions
    λ, _ = responses[(version, :rigid)]
    λFlex, _ = responses[(version, :flexible)]

    # Frequencies [Hz] of the first modes of the wing with the rigid and with the flexible tip rod
    f = sqrt.(λ[1:nModes]) / (2π)
    fFlex = sqrt.(λFlex[1:nModes]) / (2π)
    println("\nNatural frequencies of the Pazy wing S10 $version (Nastran) with the rigid and with the flexible tip rod:")
    @printf("%-4s %14s %14s %10s\n", "Mode", "Rigid [Hz]", "Flexible [Hz]", "Diff [%]")
    for i in 1:nModes
        @printf("%-4d %14.4f %14.4f %10.2f\n", i, f[i], fFlex[i], 100 * (fFlex[i] - f[i]) / f[i])
    end

    # Frequencies and mode shapes of the first modes with each tip, to be copied to AeroBeams (test/referenceData/Pazy)
    for tip in tips
        foreach(file -> println("Saved ", file), nastran_pazy_swept_modal_files(sweep, version; tip, nmodes=nModes))
    end
end
