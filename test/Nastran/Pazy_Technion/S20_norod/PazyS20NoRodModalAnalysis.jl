using BeamROM, LinearAlgebra, Printf
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

# Technion Pazy wing with 20 degrees of sweep and no tip rod modal (SOL 103) analysis with NASTRAN. The .bdf of every version are generated, and, if the results of one are not there yet, the others are still processed, so that all the .bdf to run (e.g. on a Windows VM) are listed together; run again after copying the .op2 back (S20_<version>_norod_modal.op2, in nastran/Pazy_Technion/S20/NoRod). The natural frequencies of the first modes of the wing clamped at the root are printed, and the frequencies of each are saved as S20_<version>_norod_modal_frequencies.csv (mode number, eigenvalue [(rad/s)²], frequency [Hz]) and the mode shapes as S20_<version>_norod_modal_shapes_wing.csv (translations of the sections of the spar, see extractModeShapes.py), to be copied to AeroBeams (test/referenceData/Pazy)

# Versions of the model ("LE", "TE" and/or "Baseline")
versions = ["LE", "TE", "Baseline"]

# Sweep angle [deg]: the sections (planes of constant y) of the Nastran model are shifted along x by y tan(sweep), so the axis of the spar is inclined by this angle from the y axis
sweep = 20

# Number of modes requested in the modal analysis (at least the number of modes to compare)
nd = 15
nModes = 15
nd >= nModes || error("nd = $nd must be at least $nModes")

# Eigenvalues [(rad/s)²] of the Nastran modal analysis, for each version
responses = pazy_swept_responses(version -> nastran_pazy_swept_no_rod_modal_response(sweep, version; nd), sweep, versions)

for version in versions
    λ, _ = responses[version]

    # Frequencies [Hz] of the first modes of the wing
    f = sqrt.(λ[1:nModes]) / (2π)
    println("\nNatural frequencies of the Pazy wing S20 $version without the tip rod (Nastran):")
    @printf("%-4s %14s\n", "Mode", "Frequency [Hz]")
    for i in 1:nModes
        @printf("%-4d %14.4f\n", i, f[i])
    end

    # Frequencies and mode shapes of the first modes, to be copied to AeroBeams (test/referenceData/Pazy)
    foreach(file -> println("Saved ", file), nastran_pazy_swept_no_rod_modal_files(sweep, version; nmodes=nModes))
end
