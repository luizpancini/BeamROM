using DelimitedFiles

"""
    nastran_static_response(element, mesh; L, w, t, E, ν, ρ)

Static Nastran workflow for simple beams modelled with `element` = "solid" (CHEXA), "shell" (CQUAD4) or "line" (CBAR): generates the BDF (generateBeam.py), runs Nastran (runNastran.py) and extracts the response (extractResponse.py). Files are kept in nastran/simple_beams/<element>_beam/mesh<mesh>.

Returns the (6 x 6) generalized displacements at the load node: rows = [u; θ] (translations and rotations), columns = load cases [F1, F2, F3, M1, M2, M3].

If Nastran is not available on this machine, an error is raised after listing the .bdf to run elsewhere (e.g. on a Windows VM, with `python nastran/runNastran.py <dir>`). Copy the .op2 back to that directory and call this function again. The Python interpreter is the project's `.venv`
unless the BEAMROM_PYTHON environment variable is set.
"""
function nastran_static_response(element::AbstractString, mesh::Integer; L, w, t, E, ν, ρ)
    nastranDir = @__DIR__
    repoDir = dirname(nastranDir)
    job = "$(element)_beam_static"
    caseDir = joinpath(nastranDir, "simple_beams", "$(element)_beam", "mesh$mesh")
    venvPython = joinpath(repoDir, ".venv", Sys.iswindows() ? joinpath("Scripts", "python.exe") : joinpath("bin", "python"))
    python = get(ENV, "BEAMROM_PYTHON", isfile(venvPython) ? venvPython : "python3")

    run(`$python $(joinpath(nastranDir, "generateBeam.py")) --element $element --mesh $mesh --analyses static --L $L --w $w --t $t --E $E --nu $ν --rho $ρ`)

    # Exit code 2: Nastran is not installed here
    nastranRun = run(ignorestatus(`$python $(joinpath(nastranDir, "runNastran.py")) $caseDir --jobs $job`))
    nastranRun.exitcode == 0 || error("Nastran has not produced the results (exit code $(nastranRun.exitcode)); see the messages above")

    return read_extracted_csv(`$python $(joinpath(nastranDir, "extractResponse.py")) $caseDir`)
end

# Matrix printed as CSV on the standard output by an extraction script (extractResponse.py, extractStations.py, extractSegments.py), whose messages go to the standard error. No intermediate file is written
read_extracted_csv(cmd) = readdlm(IOBuffer(read(cmd, String)), ',')

const PAZY_MODAL_JOB = "S0_modal"

function pazy_python_and_dirs()
    nastranDir = @__DIR__
    venvPython = joinpath(dirname(nastranDir), ".venv", Sys.iswindows() ? joinpath("Scripts", "python.exe") : joinpath("bin", "python"))
    python = get(ENV, "BEAMROM_PYTHON", isfile(venvPython) ? venvPython : "python3")
    return python, nastranDir, joinpath(nastranDir, "Pazy_Technion", "S0")
end

const PAZY_STRAIGHT_SWEEP = 0 # [deg] sweep angle of the straight Pazy wing models (S0 and S0_norod, unswept by definition); the swept wings (S10, S20) take theirs as an argument

"""
    pazy_reference_axis_tag(reference_axis) -> String

Label of the reference axis in the names of the Pazy files, `"RA<x>"` for the chord fraction `x`, between 0 and 1 (with "p" as decimal point, e.g. `0.43` -> `"RA0p43"`).
"""
function pazy_reference_axis_tag(reference_axis::Real)
    0 <= reference_axis <= 1 || throw(ArgumentError("reference_axis must be a fraction of the chord between 0 and 1, got $reference_axis"))
    return "RA" * replace(string(reference_axis), "." => "p")
end

"""
    pazy_mesh_tag(mesh) -> String

Label of the beam mesh in the names of the Pazy files, `"<x>mesh"`: the name of the symbol (e.g. `:nominal` -> `"nominalmesh"`, `:refined` -> `"refinedmesh"`), the number of nodes (e.g. `16` -> `"16mesh"`) or, for a vector of fractions of the span, its length and a checksum of its values (e.g. `"v5_1a2b3c4dmesh"`).
"""
pazy_mesh_tag(mesh::Symbol) = String(mesh) * "mesh"
pazy_mesh_tag(mesh::Integer) = "$(mesh)mesh"
function pazy_mesh_tag(mesh::AbstractVector{<:Real})
    checksum = UInt32(2166136261) # FNV-1a
    for byte in codeunits(join(Float64.(mesh), ","))
        checksum = (checksum ⊻ byte) * UInt32(16777619)
    end
    return "v$(length(mesh))_" * string(checksum, base=16, pad=8) * "mesh"
end

# Suffix of the names of the least squares stiffness files for the `segments` of the fit: none for :nodes (one matrix for each segment between the beam nodes, as in AeroBeams), "_global" for :global
function pazy_segments_suffix(segments::Symbol)
    segments in (:nodes, :global) || throw(ArgumentError("segments must be :nodes or :global, got $segments"))
    return segments == :nodes ? "" : "_$segments"
end

# Name of the Nastran job (static or segments) of the Technion Pazy wing for a tag: <model>_<tag>_<kind>
pazy_job(kind::AbstractString, tag::AbstractString) = "S0_$(tag)_$kind"

"""
    pazy_tip_suffix(tip) -> String

Suffix of the names of the jobs and files of the Technion Pazy wing with the tip rod rigid (`tip = :rigid`, `"_rigidtip"`) or flexible (`tip = :flexible`, `"_flextip"`). With the rigid rod, the tip section of the wing is rigid, as the stiffness matrices of the beam model require (the rod of the model is made of a rigid material, see `generatePazy.py --tip`), and with the flexible one the wing is as in the original model, for the modal analysis (the rod bends, with the natural frequencies of the wing).
"""
function pazy_tip_suffix(tip::Symbol)
    tip in (:rigid, :flexible) || throw(ArgumentError("tip must be :rigid or :flexible, got $tip"))
    return tip == :rigid ? "_rigidtip" : "_flextip"
end

"""
    pazy_tip_responses(f, tips=(:rigid, :flexible)) -> Dict

Calls `f(tip)` (e.g. `tip -> nastran_pazy_modal_response(; tip)`) for each of the tip rod `tips` (`:rigid` or `:flexible`, see `pazy_tip_suffix`) and returns the results in a dictionary indexed by tip. If it fails for one (typically because Nastran is not available here and its results are not there yet), the others are still tried, so that the .bdf files of all the tips are generated and listed together to be run elsewhere, and then an error is raised that lists the tips without results.
"""
function pazy_tip_responses(f, tips=(:rigid, :flexible))
    responses = Dict{Symbol,Any}()
    failed = Symbol[]
    for tip in tips
        try
            responses[tip] = f(tip)
        catch err
            err isa InterruptException && rethrow()
            push!(failed, tip)
            @warn "tip $tip: $(sprint(showerror, err))"
        end
    end
    isempty(failed) || error("no Nastran results for the tips $(join(failed, ", ")); see the messages above")
    return responses
end

"""
    pazy_mass_model_file(; reference_axis, mesh) -> String

Path of the CSV file with the equivalent beam mass model of the Technion Pazy wing for the given reference axis (chord fraction) and mesh (see `pazy_beam_mass_model`), which is nastran/Pazy_Technion/S0/S0_RA<x>_inertia_matrices_<mesh>.csv (see `pazy_mesh_tag`), as in AeroBeams (test/referenceData/Pazy) except for the reference axis.
"""
pazy_mass_model_file(; reference_axis::Real, mesh) = joinpath(@__DIR__, "Pazy_Technion", "S0", "S0_$(pazy_reference_axis_tag(reference_axis))_inertia_matrices_$(pazy_mesh_tag(mesh)).csv")

"""
    pazy_stiffness_file(; reference_axis, tip=:rigid) -> String

Path of the CSV file with the global stiffness matrix of the Technion Pazy wing for the given reference axis (chord fraction), which is nastran/Pazy_Technion/S0/S0_RA<x>_stiffness_matrix<suffix>.csv, with the suffix of the rigid or flexible tip rod `tip` (see `pazy_tip_suffix`). The mesh does not affect it.
"""
pazy_stiffness_file(; reference_axis::Real, tip::Symbol=:rigid) = joinpath(@__DIR__, "Pazy_Technion", "S0", "S0_$(pazy_reference_axis_tag(reference_axis))_stiffness_matrix$(pazy_tip_suffix(tip)).csv")

"""
    pazy_discretized_stiffness_file(; reference_axis, mesh, tip=:rigid) -> String

Path of the CSV file with the discretized stiffness matrices of the Technion Pazy wing for the given reference axis and mesh (see `nastran_pazy_station_response`), which is nastran/Pazy_Technion/S0/S0_RA<x>_discretized_stiffness_matrices_<mesh><suffix>.csv (see `pazy_mesh_tag` and `pazy_tip_suffix`).
"""
pazy_discretized_stiffness_file(; reference_axis::Real, mesh, tip::Symbol=:rigid) = joinpath(@__DIR__, "Pazy_Technion", "S0", "S0_$(pazy_reference_axis_tag(reference_axis))_discretized_stiffness_matrices_$(pazy_mesh_tag(mesh))$(pazy_tip_suffix(tip)).csv")

"""
    pazy_clamped_segment_stiffness_file(; reference_axis, mesh, intervals, method, tip=:rigid) -> String

Path of the CSV file with the clamped segment stiffness matrices of the Technion Pazy wing for the given reference axis and mesh, which is nastran/Pazy_Technion/S0/S0_RA<x>_clamped_segment_stiffness_matrices_<mesh>_intervals<intervals>_<method><suffix>.csv (see `pazy_mesh_tag` and `pazy_tip_suffix`).
"""
pazy_clamped_segment_stiffness_file(; reference_axis::Real, mesh, intervals::Integer, method, tip::Symbol=:rigid) = joinpath(@__DIR__, "Pazy_Technion", "S0", "S0_$(pazy_reference_axis_tag(reference_axis))_clamped_segment_stiffness_matrices_$(pazy_mesh_tag(mesh))_intervals$(intervals)_$method$(pazy_tip_suffix(tip)).csv")

"""
    pazy_least_squares_stiffness_file(; reference_axis, mesh, segments, tip=:rigid) -> String

Path of the CSV file with the stiffness matrices of the Technion Pazy wing fitted by `least_squares_stiffness_matrices` (`segments` = `:nodes` for one for each segment between the beam nodes, `:global` for a single one) for the given reference axis and mesh, which is nastran/Pazy_Technion/S0/S0_RA<x>_stiffness_matrices_<mesh><suffix>.csv (see `pazy_mesh_tag` and `pazy_tip_suffix`), as in AeroBeams (test/referenceData/Pazy) but with the reference axis. The file of `:global` has the suffix `_global` (see `pazy_segments_suffix`).
"""
pazy_least_squares_stiffness_file(; reference_axis::Real, mesh, segments::Symbol, tip::Symbol=:rigid) = joinpath(@__DIR__, "Pazy_Technion", "S0", "S0_$(pazy_reference_axis_tag(reference_axis))_stiffness_matrices_$(pazy_mesh_tag(mesh))$(pazy_tip_suffix(tip))$(pazy_segments_suffix(segments)).csv")

# Exit code 2: Nastran is not installed here
function run_pazy_nastran(python, nastranDir, caseDir, job)
    nastranRun = run(ignorestatus(`$python $(joinpath(nastranDir, "runNastran.py")) $caseDir --jobs $job`))
    nastranRun.exitcode == 0 || error("Nastran has not produced the results (exit code $(nastranRun.exitcode)); see the messages above")
end

"""
    nastran_pazy_static_response(; reference_axis, tip=:rigid)

Static Nastran workflow for the Technion Pazy wing (nastran/Pazy_Technion/S0/S0.bdf), whose beam is the aluminum spar: generates the static BDF (generatePazy.py), runs Nastran (runNastran.py) and extracts the response (extractResponse.py). Files are kept in nastran/Pazy_Technion/S0. The load node is at the end of the beam axis, at `reference_axis` (fraction of the chord from the minimum x of the shell grids, between 0 and 1, as in `pazy_beam_mass_model`) of the spar. Each reference axis has its own files, named S0_RA<x>_static<suffix>.* (see `pazy_reference_axis_tag`), so the results of the others are kept. `tip` is `:rigid` (default, the tip rod is rigid, so the tip section of the wing is, as the stiffness matrices require) or `:flexible` (the rod of the model, which bends), and its suffix (`pazy_tip_suffix`) is in the names, so that the results of both are kept.

Returns `(T, L)`: `T` is the (6 x 6) generalized displacements at the load node, in the beam frame (x1 along the span from the root to the tip, x2 = -x of the Nastran model, x3 = z), rows = [u; θ] and columns = load cases [F1, F2, F3, M1, M2, M3], and `L` is the length of the beam [m] (from the root section to the end of the spar).

If Nastran is not available on this machine, an error is raised after listing the .bdf to run elsewhere (e.g. on a Windows VM, with `python nastran/runNastran.py "nastran/Pazy_Technion/S0" --jobs S0_RA0p43_static_rigidtip`). Copy the .op2 file back to that directory and call this function again. The Python interpreter is the project's `.venv` unless the BEAMROM_PYTHON environment variable is set.
"""
function nastran_pazy_static_response(; reference_axis::Real, tip::Symbol=:rigid)
    python, nastranDir, caseDir = pazy_python_and_dirs()
    tag = pazy_reference_axis_tag(reference_axis)
    job = pazy_job("static", tag) * pazy_tip_suffix(tip)

    run(`$python $(joinpath(nastranDir, "generatePazy.py")) --sweep $PAZY_STRAIGHT_SWEEP --analyses static --tag $tag --chord-fraction $reference_axis --tip $tip`)
    run_pazy_nastran(python, nastranDir, caseDir, job)

    L = parse(Float64, match(r"\"L\":\s*([-+0-9.eE]+)", read(joinpath(caseDir, "$job.meta.json"), String)).captures[1])
    T = read_extracted_csv(`$python $(joinpath(nastranDir, "extractResponse.py")) $caseDir --job $job`)

    return T, L
end

"""
    nastran_pazy_station_response(y; reference_axis, tol=1e-5, twist_center=:auto, tip=:rigid)

Static Nastran workflow for the Technion Pazy wing (nastran/Pazy_Technion/S0/S0.bdf) to get the generalized displacements of the nodes of a discretized beam, whose axis is the one of `nastran_pazy_static_response` (the loads are the same: the six unit loads at the end of the axis, with the root section clamped). `y` are the spanwise positions (y of the Nastran model) of the nodes, from the root to the end of the spar, e.g. `[x[2] for x in spar_reference_axis(read_gfem(bdf), :nominal; reference_axis)]`, which must be at rows of grids of the spar (within `tol`).

The generalized displacements of a node are the rigid body motion of the section, fitted to the displacements of all the spar grids at that y, referred to the axis (extractStations.py). It uses the job of `nastran_pazy_static_response`, so no extra Nastran run is needed, and no file is written: the script prints the response (see `read_extracted_csv`).

The translations and the rotations about the span and normal axes come from the translations of all the spar grids at the node, and the rotation about the chord axis from the rotations of the grids about it (the normal of the plate), which vary linearly along the chord with the twist rate (warping of the section). The rigid body rotation of the section is the value of that fit at the twist center, where the warping vanishes, and `twist_center` sets that position: `:auto` (default) finds it from the response (the position at which the rotation about the normal axis under the moment about the span equals the rotation about the span under the moment about the normal axis, as the reciprocity of the beam model requires, which is the same for every reference axis; it is printed), a number is a fraction of the chord (as `reference_axis`) and `:axis` takes it at the reference axis. With `:axis` the response, and hence the stiffness matrices, depend on the position of the axis through the warping of the section, which the beam model cannot represent (e.g. its shear center shifts towards the axis); with the twist center they do not (the beam is the same one for any reference axis, with the matrices at another point related by a rigid transformation).

Returns `(T, L)`: `T` is a vector with the (6 x 6) generalized displacements of each node, in the beam frame (x1 along the span, x2 = -x of the Nastran model, x3 = z), rows = [u; θ] and columns = load cases [F1, F2, F3, M1, M2, M3] (zeros at the clamped root), and `L` is the length of the beam [m] (from the root section to the end of the spar).

If Nastran is not available on this machine, an error is raised after listing the .bdf to run elsewhere (see `nastran_pazy_static_response`). The Python interpreter is the project's `.venv` unless the BEAMROM_PYTHON environment variable is set.
"""
function nastran_pazy_station_response(y::AbstractVector{<:Real}; reference_axis::Real, tol::Real=1e-5, twist_center::Union{Symbol,Real}=:auto, tip::Symbol=:rigid)
    twist_center isa Real || twist_center in (:auto, :axis) || throw(ArgumentError("twist_center must be :auto, :axis or a fraction of the chord, got $twist_center"))
    python, nastranDir, caseDir = pazy_python_and_dirs()
    axisTag = pazy_reference_axis_tag(reference_axis)
    job = pazy_job("static", axisTag) * pazy_tip_suffix(tip)
    generate = joinpath(nastranDir, "generatePazy.py")
    extract = joinpath(nastranDir, "extractStations.py")

    run(`$python $generate --sweep $PAZY_STRAIGHT_SWEEP --analyses static --tag $axisTag --chord-fraction $reference_axis --tip $tip`)
    run_pazy_nastran(python, nastranDir, caseDir, job)

    L = parse(Float64, match(r"\"L\":\s*([-+0-9.eE]+)", read(joinpath(caseDir, "$job.meta.json"), String)).captures[1])
    Tall = read_extracted_csv(`$python $extract $caseDir --job $job --stations $y --tol $tol --twist-center $twist_center`)
    T = [Tall[6*(k-1)+1:6*k, :] for k in eachindex(y)]

    return T, L
end

"""
    nastran_pazy_segment_response(y; intervals=1, reference_axis, mesh, tol=1e-5, tip=:rigid)

Static Nastran workflow for the Technion Pazy wing (nastran/Pazy_Technion/S0/S0.bdf) to get the generalized displacements of windows of `intervals` consecutive segments of a discretized beam analyzed as cantilevers: in a single job (S0_RA<x>_<mesh>_segments<suffix>, with the suffix of the rigid or flexible tip rod `tip`, see `nastran_pazy_static_response`; generatePazy.py --analyses segments, extractSegments.py) with six load cases for each window `w` (`w = 1, ..., n - intervals`), the section of the wing at the node `w` is clamped and the six unit loads are applied at the node `w + intervals`, whose displacements are taken. The rest of the wing is kept. The same BDF is overwritten for each `intervals`; when the content of the BDF changes (`intervals`, `y`, `reference_axis`) the .op2 and the other outputs of the previous job are deleted, and Nastran has to be run again. `y` are the spanwise positions (y of the Nastran model) of the nodes, from the root to the end of the spar, as in `nastran_pazy_station_response`, and `mesh` is the one that gives them (see `pazy_mesh_tag`), which labels the job: each combination of reference axis and mesh has its own. The nodes are the grids on the beam axis tied by an RBE3 to the spar grids at that y. The loads are applied where the displacements are taken, so the displacement-force matrices of the segments are reciprocal.

Returns `(T, L)`: `T` is a vector with the `length(y) - intervals` (6 x 6) generalized displacements of the downstream node of each window, in the beam frame (x1 along the span, x2 = -x of the Nastran model, x3 = z), rows = [u; θ] and columns = load cases [F1, F2, F3, M1, M2, M3], and `L` is the length of the beam [m] (from the root section to the end of the spar).

If Nastran is not available on this machine, an error is raised after listing the .bdf to run elsewhere (e.g. on a Windows VM, with `python nastran/runNastran.py "nastran/Pazy_Technion/S0" --jobs S0_RA0p43_nominalmesh_segments`). Copy the .op2 file back to that directory and call this function again. The Python interpreter is the project's `.venv` unless the BEAMROM_PYTHON environment variable is set.
"""
function nastran_pazy_segment_response(y::AbstractVector{<:Real}; intervals::Integer=1, reference_axis::Real, mesh, tol::Real=1e-5, tip::Symbol=:rigid)
    python, nastranDir, caseDir = pazy_python_and_dirs()
    tag = "$(pazy_reference_axis_tag(reference_axis))_$(pazy_mesh_tag(mesh))"
    job = pazy_job("segments", tag) * pazy_tip_suffix(tip)

    run(`$python $(joinpath(nastranDir, "generatePazy.py")) --sweep $PAZY_STRAIGHT_SWEEP --analyses segments --tag $tag --chord-fraction $reference_axis --intervals $intervals --stations $y --tol $tol --tip $tip`)
    run_pazy_nastran(python, nastranDir, caseDir, job)

    L = parse(Float64, match(r"\"L\":\s*([-+0-9.eE]+)", read(joinpath(caseDir, "$job.json"), String)).captures[1])
    Tall = read_extracted_csv(`$python $(joinpath(nastranDir, "extractSegments.py")) $caseDir --job $job`)
    T = [Tall[6*(k-1)+1:6*k, :] for k in 1:length(y)-intervals]

    return T, L
end

"""
    nastran_pazy_modal_response(; nd=10, tip=:flexible)

Modal Nastran workflow for the Technion Pazy wing (nastran/Pazy_Technion/S0/S0.bdf): generates the modal BDF (generatePazy.py), runs Nastran (runNastran.py) and extracts the eigenvalues (extractModes.py). Files are kept in nastran/Pazy_Technion/S0. `nd` is the number of modes requested. The beam reference axis does not affect the modal analysis, so its files have no reference axis label, but they have the suffix of the rigid or flexible tip rod `tip` (S0_modal<suffix>.*, see `nastran_pazy_static_response`; the default is `:flexible`, the wing as it is, and the `:rigid` one has a tip like that of the beam model for the stiffness).

Returns `(λ, L)`: the `nd` eigenvalues [(rad/s)²] of the wing clamped at the root, in ascending order (frequencies in Hz are `sqrt.(λ) / 2π`), and the length of the beam [m] (from the root section to the end of the spar).

If Nastran is not available on this machine, an error is raised after listing the .bdf to run elsewhere (e.g. on a Windows VM, with `python nastran/runNastran.py "nastran/Pazy_Technion/S0" --jobs S0_modal_flextip`). Copy the .op2 file back to that directory and call this function again. The Python interpreter is the project's `.venv` unless the BEAMROM_PYTHON environment variable is set.
"""
function nastran_pazy_modal_response(; nd::Integer=10, tip::Symbol=:flexible)
    python, nastranDir, caseDir = pazy_python_and_dirs()
    job = PAZY_MODAL_JOB * pazy_tip_suffix(tip)
    frequencies = "S0_modal_frequencies$(pazy_tip_suffix(tip)).csv"

    run(`$python $(joinpath(nastranDir, "generatePazy.py")) --sweep $PAZY_STRAIGHT_SWEEP --analyses modal --nd $nd --tip $tip`)
    run_pazy_nastran(python, nastranDir, caseDir, job)
    run(`$python $(joinpath(nastranDir, "extractModes.py")) $caseDir --job $job --out $frequencies`)

    # Length of the beam, written by generatePazy.py next to the modal deck
    L = parse(Float64, match(r"\"L\":\s*([-+0-9.eE]+)", read(joinpath(caseDir, "$job.json"), String)).captures[1])

    return readdlm(joinpath(caseDir, frequencies), ',')[:, 2], L
end

"""
    nastran_pazy_modal_shapes(; tip=:flexible, nmodes=15) -> Vector{String}

Extracts the mode shapes of the first `nmodes` modes of the modal analysis of the Technion Pazy wing with the tip rod `tip` (the .op2 of `nastran_pazy_modal_response`, S0_modal<suffix>.op2) with extractModeShapes.py, to compare them with those of the beam model (modal assurance criterion): the translations of the leading and trailing edge points of the sections of the spar and of the grids of the tip rod. Saves them as S0_modal_shapes<suffix>_wing.csv and S0_modal_shapes<suffix>_rod.csv in nastran/Pazy_Technion/S0 (to be copied to AeroBeams, test/referenceData/Pazy) and returns their paths.
"""
function nastran_pazy_modal_shapes(; tip::Symbol=:flexible, nmodes::Integer=15)
    python, nastranDir, caseDir = pazy_python_and_dirs()
    suffix = pazy_tip_suffix(tip)
    job = PAZY_MODAL_JOB * suffix
    name = "S0_modal_shapes$suffix"

    run(`$python $(joinpath(nastranDir, "extractModeShapes.py")) $caseDir --job $job --nmodes $nmodes --name $name`)

    return [joinpath(caseDir, "$(name)_$kind.csv") for kind in ("wing", "rod")]
end

# ------------------------------------------------------------------------------------------------------------
# Swept Technion Pazy wings (S10 and S20, with 10 and 20 degrees of sweep; versions LE, TE and Baseline). The sweep angle is never assumed: every function takes it, as its first argument
# ------------------------------------------------------------------------------------------------------------

const PAZY_SWEPT_VERSIONS = ("LE", "TE", "Baseline")

"""
    pazy_swept_label(sweep) -> String

Name of the swept Technion Pazy wing with `sweep` degrees of sweep (an integer number of degrees), "S<sweep>" (e.g. "S10"), which names its folder in nastran/Pazy_Technion and its files.
"""
function pazy_swept_label(sweep::Real)
    isinteger(sweep) || throw(ArgumentError("the sweep angle of a swept Pazy wing must be an integer number of degrees, got $sweep"))
    return "S$(Int(sweep))"
end

# Folder of a version of a swept wing, nastran/Pazy_Technion/S<sweep>/<version>
function pazy_swept_dir(sweep::Real, version::AbstractString)
    version in PAZY_SWEPT_VERSIONS || throw(ArgumentError("unknown version $version (use $(join(PAZY_SWEPT_VERSIONS, ", ", " or ")))"))
    return joinpath(@__DIR__, "Pazy_Technion", pazy_swept_label(sweep), version)
end

# Name of the prepared .bdf file and of the Nastran jobs of a version of a swept wing: S<sweep>_LE, S<sweep>_TE and S<sweep>_baseline
pazy_swept_name(sweep::Real, version::AbstractString) = "$(pazy_swept_label(sweep))_$(version == "Baseline" ? "baseline" : version)"

"""
    pazy_swept_mass_model_file(sweep, version; reference_axis, mesh) -> String

Path of the CSV file with the equivalent beam mass model of the swept wing with `sweep` degrees of sweep (S10, S20), `version`, for the given reference axis and mesh (see `pazy_beam_mass_model`), which is nastran/Pazy_Technion/S<sweep>/<version>/S<sweep>_<version>_RA<x>_inertia_matrices_<mesh>.csv (see `pazy_mesh_tag`; the version is as in `pazy_swept_name`, e.g. S10_baseline), as in AeroBeams (test/referenceData/Pazy) except for the reference axis.
"""
pazy_swept_mass_model_file(sweep::Real, version::AbstractString; reference_axis::Real, mesh) = joinpath(pazy_swept_dir(sweep, version), "$(pazy_swept_name(sweep, version))_$(pazy_reference_axis_tag(reference_axis))_inertia_matrices_$(pazy_mesh_tag(mesh)).csv")

# Elements of the tip rod of the swept models that a model without the rod leaves out: the CBEAM properties of the rod (13) and of the breaking beam (24, LE and TE), and the CONM2 of the fixed adapter and of the FPS weights (the calibration mass is part of the wing and stays). Their ids depend on the wing: S10 has the fixed adapter 205 and the FPS weights 433 (LE and TE); S20 has the FPS weights 203 (LE) and 433 (TE), and the fixed adapter 204 (LE and TE), and its 205 is the calibration mass
const PAZY_SWEPT_ROD_BEAM_PIDS = (13, 24)
const PAZY_SWEPT_ROD_CONM2 = Dict(10 => (205, 433), 20 => (203, 204, 433))

"""
    pazy_swept_rod_conm2(sweep) -> Tuple

Ids of the CONM2 of the fixed adapter and of the FPS weights of the swept wing with `sweep` degrees of sweep, which a model without the tip rod leaves out (see `PAZY_SWEPT_ROD_CONM2`).
"""
function pazy_swept_rod_conm2(sweep::Real)
    haskey(PAZY_SWEPT_ROD_CONM2, sweep) || throw(ArgumentError("the CONM2 of the tip rod of the swept wing with $sweep degrees of sweep are not known (see PAZY_SWEPT_ROD_CONM2)"))
    return PAZY_SWEPT_ROD_CONM2[sweep]
end

# Folder of the decks and of the Nastran files of a swept wing without the tip rod, nastran/Pazy_Technion/S<sweep>/NoRod (the decks are S<sweep>_<version>_norod.bdf)
const PAZY_SWEPT_NO_ROD_DIR = "NoRod"

pazy_swept_no_rod_dir(sweep::Real, version::AbstractString) = joinpath(dirname(pazy_swept_dir(sweep, version)), PAZY_SWEPT_NO_ROD_DIR)

# Name of the deck and of the Nastran jobs of a version of a swept wing without the tip rod: S<sweep>_LE_norod, S<sweep>_TE_norod and S<sweep>_baseline_norod
pazy_swept_no_rod_name(sweep::Real, version::AbstractString) = pazy_swept_name(sweep, version) * "_norod"

pazy_swept_no_rod_bdf(sweep::Real, version::AbstractString) = joinpath(pazy_swept_no_rod_dir(sweep, version), "$(pazy_swept_no_rod_name(sweep, version)).bdf")

function pazy_swept_no_rod_dirs(sweep::Real, version::AbstractString)
    python, nastranDir, _ = pazy_python_and_dirs()
    return python, nastranDir, pazy_swept_no_rod_dir(sweep, version), pazy_swept_no_rod_name(sweep, version)
end

"""
    pazy_swept_no_rod_mass_model_file(sweep, version; reference_axis, mesh) -> String

Path of the CSV file with the equivalent beam mass model of the swept wing with `sweep` degrees of sweep, `version`, without the tip rod (nor the adapter and the FPS weights, see `PAZY_SWEPT_ROD_BEAM_PIDS` and `pazy_swept_rod_conm2`; the calibration mass is part of the wing and stays), for the given reference axis and mesh, which is nastran/Pazy_Technion/S<sweep>/NoRod/S<sweep>_<version>_norod_RA<x>_inertia_matrices_<mesh>.csv (see `pazy_mesh_tag`; the version is as in `pazy_swept_name`).
"""
pazy_swept_no_rod_mass_model_file(sweep::Real, version::AbstractString; reference_axis::Real, mesh) = joinpath(pazy_swept_no_rod_dir(sweep, version), "$(pazy_swept_no_rod_name(sweep, version))_$(pazy_reference_axis_tag(reference_axis))_inertia_matrices_$(pazy_mesh_tag(mesh)).csv")

"""
    pazy_swept_stiffness_file(sweep, version; reference_axis, tip=:rigid) -> String

Path of the CSV file with the global stiffness matrix of the swept wing with `sweep` degrees of sweep, `version`, for the given reference axis (chord fraction or `nothing`), which is nastran/Pazy_Technion/S<sweep>/<version>/S<sweep>_<version>_RA<x>_stiffness_matrix<suffix>.csv, with the suffix of the rigid or flexible tip rod `tip` (see `pazy_tip_suffix`). The mesh does not affect it.
"""
pazy_swept_stiffness_file(sweep::Real, version::AbstractString; reference_axis::Real, tip::Symbol=:rigid) = joinpath(pazy_swept_dir(sweep, version), "$(pazy_swept_label(sweep))_$(version)_$(pazy_reference_axis_tag(reference_axis))_stiffness_matrix$(pazy_tip_suffix(tip)).csv")

"""
    pazy_swept_least_squares_stiffness_file(sweep, version; reference_axis, mesh, segments, tip=:rigid) -> String

Path of the CSV file with the stiffness matrices of the swept wing with `sweep` degrees of sweep, `version`, fitted by `least_squares_stiffness_matrices` (`segments` = `:nodes` or `:global`, see `pazy_least_squares_stiffness_file`) for the given reference axis and mesh, which is nastran/Pazy_Technion/S<sweep>/<version>/S<sweep>_<version>_RA<x>_stiffness_matrices_<mesh><suffix>.csv (see `pazy_mesh_tag`; the version is as in `pazy_swept_name`, e.g. S10_baseline), with the suffix of the rigid or flexible tip rod `tip` (see `pazy_tip_suffix`), and `_global` after it for `segments = :global` (see `pazy_segments_suffix`).
"""
pazy_swept_least_squares_stiffness_file(sweep::Real, version::AbstractString; reference_axis::Real, mesh, segments::Symbol, tip::Symbol=:rigid) = joinpath(pazy_swept_dir(sweep, version), "$(pazy_swept_name(sweep, version))_$(pazy_reference_axis_tag(reference_axis))_stiffness_matrices_$(pazy_mesh_tag(mesh))$(pazy_tip_suffix(tip))$(pazy_segments_suffix(segments)).csv")

"""
    pazy_swept_no_rod_stiffness_file(sweep, version; reference_axis) -> String

Path of the CSV file with the global stiffness matrix of the swept wing with `sweep` degrees of sweep, `version`, without the tip rod, for the given reference axis (chord fraction), which is nastran/Pazy_Technion/S<sweep>/NoRod/S<sweep>_<version>_norod_RA<x>_stiffness_matrix.csv.
"""
pazy_swept_no_rod_stiffness_file(sweep::Real, version::AbstractString; reference_axis::Real) = joinpath(pazy_swept_no_rod_dir(sweep, version), "$(pazy_swept_label(sweep))_$(version)_norod_$(pazy_reference_axis_tag(reference_axis))_stiffness_matrix.csv")

"""
    pazy_swept_no_rod_least_squares_stiffness_file(sweep, version; reference_axis, mesh, segments) -> String

Path of the CSV file with the stiffness matrices of the swept wing with `sweep` degrees of sweep, `version`, without the tip rod, fitted by `least_squares_stiffness_matrices` (`segments` = `:nodes` or `:global`) for the given reference axis and mesh, which is nastran/Pazy_Technion/S<sweep>/NoRod/S<sweep>_<version>_norod_RA<x>_stiffness_matrices_<mesh>.csv (see `pazy_mesh_tag`; the version is as in `pazy_swept_name`), and `_global` after it for `segments = :global` (see `pazy_segments_suffix`).
"""
pazy_swept_no_rod_least_squares_stiffness_file(sweep::Real, version::AbstractString; reference_axis::Real, mesh, segments::Symbol) = joinpath(pazy_swept_no_rod_dir(sweep, version), "$(pazy_swept_no_rod_name(sweep, version))_$(pazy_reference_axis_tag(reference_axis))_stiffness_matrices_$(pazy_mesh_tag(mesh))$(pazy_segments_suffix(segments)).csv")


function pazy_swept_dirs(sweep::Real, version::AbstractString)
    python, nastranDir, _ = pazy_python_and_dirs()
    return python, nastranDir, pazy_swept_dir(sweep, version), pazy_swept_name(sweep, version)
end

"""
    prepare_pazy_swept(sweep, version)

Bulk data of the swept Pazy wing with `sweep` degrees of sweep (S10, S20), `version` = "LE", "TE" or "Baseline" (the position of the tip weights, nastran/Pazy_Technion/S<sweep>/<version>/S<sweep>_<version>.nas, which are not modified; the S10 Baseline deck is nastran/Pazy_Technion/S10/Baseline/S10_baseline.bdf, and the S20 one is made from the LE deck): runs prepareSweptPazy.py, which removes the aerodynamic and flutter data of the Femap deck, leaves the materials as they are (flexible tip rod and breaking beam: the rigid nylon of the tip is set in the decks of the analyses, see `nastran_pazy_swept_static_response`), ties the roll of the leading and trailing edge stringers to the spar at the ribs, and writes the prepared .bdf (S<sweep>_<version>.bdf; for the S10 Baseline, its own S10_baseline.bdf is prepared in place) only if its content changes, so that the decks and the Nastran results that come from it are kept. Returns the path of the .bdf.
"""
function prepare_pazy_swept(sweep::Real, version::AbstractString)
    python, nastranDir, caseDir, name = pazy_swept_dirs(sweep, version)
    run(`$python $(joinpath(nastranDir, "prepareSweptPazy.py")) --sweeps $(Int(sweep)) --versions $version`)
    return joinpath(caseDir, "$name.bdf")
end

"""
    pazy_swept_responses(f, sweep, versions) -> Dict

Calls `f(version)` (e.g. `version -> nastran_pazy_swept_static_response(sweep, version; reference_axis)`) for each of the `versions` of the swept wing with `sweep` degrees of sweep and returns the results in a dictionary indexed by version. If it fails for one (typically because Nastran is not available here and its results are not there yet), the others are still tried, so that the .bdf files of all the versions are generated and listed together to be run elsewhere, and then an error is raised that lists the versions without results.
"""
function pazy_swept_responses(f, sweep::Real, versions)
    responses = Dict{String,Any}()
    failed = String[]
    for version in versions
        try
            responses[version] = f(version)
        catch err
            err isa InterruptException && rethrow()
            push!(failed, version)
            @warn "$(pazy_swept_label(sweep)) $version: $(sprint(showerror, err))"
        end
    end
    isempty(failed) || error("no Nastran results for the versions $(join(failed, ", ")); see the messages above")
    return responses
end

"""
    pazy_swept_tip_responses(f, sweep, versions, tips=(:rigid, :flexible)) -> Dict

Calls `f(version, tip)` (e.g. `(version, tip) -> nastran_pazy_swept_modal_response(sweep, version; tip)`) for each of the `versions` of the swept wing with `sweep` degrees of sweep and each of the tip rod `tips` (`:rigid` or `:flexible`, see `pazy_tip_suffix`) and returns the results in a dictionary indexed by `(version, tip)`. If it fails for one (typically because Nastran is not available here and its results are not there yet), the others are still tried, so that the .bdf files of all of them are generated and listed together to be run elsewhere, and then an error is raised that lists those without results.
"""
function pazy_swept_tip_responses(f, sweep::Real, versions, tips=(:rigid, :flexible))
    responses = Dict{Tuple{String,Symbol},Any}()
    failed = String[]
    for version in versions, tip in tips
        try
            responses[(version, tip)] = f(version, tip)
        catch err
            err isa InterruptException && rethrow()
            push!(failed, "$version ($tip)")
            @warn "$(pazy_swept_label(sweep)) $version, tip $tip: $(sprint(showerror, err))"
        end
    end
    isempty(failed) || error("no Nastran results for $(join(failed, ", ")); see the messages above")
    return responses
end

"""
    nastran_pazy_swept_static_response(sweep, version; reference_axis, tip=:rigid)

Static Nastran workflow for the swept Technion Pazy wing with `sweep` degrees of sweep (S10, S20), `version` = "LE", "TE" or "Baseline" (see `prepare_pazy_swept`; files are kept in nastran/Pazy_Technion/S<sweep>/<version>): generates the static BDF (generatePazy.py), runs Nastran (runNastran.py) and extracts the response (extractResponse.py). As in `nastran_pazy_static_response`, the wing is clamped at the root section and the six unit loads are applied at the load node, at the end of the beam axis, at `reference_axis` (fraction of the chord between 0 and 1, as in `pazy_beam_mass_model`) of the spar. The axis is inclined by `sweep` [deg] from the y axis of the Nastran model (the sections are shifted along x by `y tan(sweep)`), at a constant fraction of the chord along the span. Each reference axis has its own files, named S<sweep>_<version>_RA<x>_static<suffix>.* (see `pazy_reference_axis_tag`), so the results of the others are kept. `tip` is `:rigid` (default: the tip rod, property 13, and the breaking beam, property 24, of LE and TE, have a rigid nylon, E = 1e15 Pa with the Poisson's ratio and density of the Nylon 12, so the tip section of the wing is, as the stiffness matrices require) or `:flexible` (the materials of the model, so the rod bends), and its suffix (`pazy_tip_suffix`) is in the names, so that the results of both are kept.

Returns `(T, L)`: `T` is the (6 x 6) generalized displacements at the load node, in the beam frame rotated by `sweep` about z (x1 along the axis, from the root to the end of the spar, x2 = (-cos(sweep), sin(sweep), 0) and x3 = z of the Nastran model), rows = [u; θ] and columns = load cases [F1, F2, F3, M1, M2, M3], and `L` is the length of the beam axis [m] (from the root section to the end of the spar).

If Nastran is not available on this machine, an error is raised after listing the .bdf to run elsewhere (e.g. on a Windows VM, with `python nastran/runNastran.py "nastran/Pazy_Technion/S10/Baseline" --jobs S10_baseline_RA0p43_static_rigidtip`). Copy the .op2 file back to that directory and call this function again. The Python interpreter is the project's `.venv` unless the BEAMROM_PYTHON environment variable is set.
"""
function nastran_pazy_swept_static_response(sweep::Real, version::AbstractString; reference_axis::Real, tip::Symbol=:rigid)
    python, nastranDir, caseDir, name = pazy_swept_dirs(sweep, version)
    bdf = prepare_pazy_swept(sweep, version)
    tag = pazy_reference_axis_tag(reference_axis)
    job = "$(name)_$(tag)_static$(pazy_tip_suffix(tip))"

    run(`$python $(joinpath(nastranDir, "generatePazy.py")) --bdf $bdf --analyses static --sweep $sweep --tag $tag --chord-fraction $reference_axis --tip $tip --rigid-pids $(collect(PAZY_SWEPT_ROD_BEAM_PIDS))`)
    run_pazy_nastran(python, nastranDir, caseDir, job)

    L = parse(Float64, match(r"\"L\":\s*([-+0-9.eE]+)", read(joinpath(caseDir, "$job.meta.json"), String)).captures[1])
    T = read_extracted_csv(`$python $(joinpath(nastranDir, "extractResponse.py")) $caseDir --job $job`)

    return T, L
end

"""
    nastran_pazy_swept_station_response(sweep, version, y; reference_axis, tol=1e-5, twist_center=:auto, tip=:rigid)

Generalized displacements of the nodes of a discretized beam of the swept Technion Pazy wing with `sweep` degrees of sweep (S10, S20), `version` = "LE", "TE" or "Baseline", under the six unit loads at the end of the beam axis, from the job of `nastran_pazy_swept_static_response` (which is generated and run if needed, so no extra Nastran run is required): the rigid body motion of each section, fitted to the displacements of all the spar grids at that y and referred to the point of the axis there (extractStations.py, as in `nastran_pazy_station_response`). `y` are the spanwise positions (y of the Nastran model) of the nodes, from the root to the end of the spar, e.g. `[x[2] for x in spar_reference_axis(read_gfem(bdf), :nominal; reference_axis, sweep)]`, which must be at rows of grids of the spar (within `tol`). The sections are those of constant y, so the nodes are at the distances `(y .- y[1]) ./ cosd(sweep)` from the root along the inclined axis. No file is written: the script prints the response (see `read_extracted_csv`). The rotation about the chord axis of the sections is taken at the twist center, as in `nastran_pazy_station_response` (`twist_center` = `:auto`, the default, a fraction of the chord or `:axis`), so the response does not depend on `reference_axis`: the twist center of the swept wings is found the same way (the position at which the rotations about the normal and span axes under the unit moments are reciprocal), which for the S10 Baseline is at 0.449 of the chord for every reference axis.

Returns `(T, L)`: `T` is a vector with the (6 x 6) generalized displacements of each node, in the beam frame rotated by `sweep` about z (as in `nastran_pazy_swept_static_response`), rows = [u; θ] and columns = load cases [F1, F2, F3, M1, M2, M3] (zeros at the clamped root), and `L` is the length of the beam axis [m].
"""
function nastran_pazy_swept_station_response(sweep::Real, version::AbstractString, y::AbstractVector{<:Real}; reference_axis::Real, tol::Real=1e-5, twist_center::Union{Symbol,Real}=:auto, tip::Symbol=:rigid)
    twist_center isa Real || twist_center in (:auto, :axis) || throw(ArgumentError("twist_center must be :auto, :axis or a fraction of the chord, got $twist_center"))
    python, nastranDir, caseDir, name = pazy_swept_dirs(sweep, version)
    _, L = nastran_pazy_swept_static_response(sweep, version; reference_axis, tip)
    job = "$(name)_$(pazy_reference_axis_tag(reference_axis))_static$(pazy_tip_suffix(tip))"

    Tall = read_extracted_csv(`$python $(joinpath(nastranDir, "extractStations.py")) $caseDir --job $job --stations $y --tol $tol --twist-center $twist_center`)
    T = [Tall[6*(k-1)+1:6*k, :] for k in eachindex(y)]

    return T, L
end

"""
    nastran_pazy_swept_modal_response(sweep, version; nd=10, tip=:flexible)

Modal Nastran workflow for the swept Technion Pazy wing with `sweep` degrees of sweep (S10, S20), `version` = "LE", "TE" or "Baseline" (see `prepare_pazy_swept`; files are kept in nastran/Pazy_Technion/S<sweep>/<version>): generates the modal BDF (generatePazy.py), runs Nastran (runNastran.py) and extracts the eigenvalues (extractModes.py). `nd` is the number of modes requested. The beam reference axis does not affect the modal analysis, `sweep` [deg] only sets the length of the beam axis that is returned. `tip` is `:flexible` (default, the wing as it is, with the rod and the breaking beam of the model) or `:rigid` (rigid nylon in the properties 13 and 24, see `nastran_pazy_swept_static_response`), and its suffix (`pazy_tip_suffix`) is in the names of the files (S<sweep>_<version>_modal<suffix>.*), so that the results of both are kept.

Returns `(λ, L)`: the `nd` eigenvalues [(rad/s)²] of the wing clamped at the root, in ascending order (frequencies in Hz are `sqrt.(λ) / 2π`), and the length of the beam axis [m] (from the root section to the end of the spar).

If Nastran is not available on this machine, an error is raised after listing the .bdf to run elsewhere (e.g. on a Windows VM, with `python nastran/runNastran.py "nastran/Pazy_Technion/S10/LE" --jobs S10_LE_modal_flextip`). Copy the .op2 file back to that directory and call this function again. The Python interpreter is the project's `.venv` unless the BEAMROM_PYTHON environment variable is set.
"""
function nastran_pazy_swept_modal_response(sweep::Real, version::AbstractString; nd::Integer=10, tip::Symbol=:flexible)
    python, nastranDir, caseDir, name = pazy_swept_dirs(sweep, version)
    bdf = prepare_pazy_swept(sweep, version)
    job = "$(name)_modal$(pazy_tip_suffix(tip))"
    frequencies = "$(name)_modal_frequencies$(pazy_tip_suffix(tip)).csv"

    run(`$python $(joinpath(nastranDir, "generatePazy.py")) --bdf $bdf --analyses modal --nd $nd --sweep $sweep --tip $tip --rigid-pids $(collect(PAZY_SWEPT_ROD_BEAM_PIDS))`)
    run_pazy_nastran(python, nastranDir, caseDir, job)
    run(`$python $(joinpath(nastranDir, "extractModes.py")) $caseDir --job $job --out $frequencies`)

    L = parse(Float64, match(r"\"L\":\s*([-+0-9.eE]+)", read(joinpath(caseDir, "$job.json"), String)).captures[1])

    return readdlm(joinpath(caseDir, frequencies), ',')[:, 2], L
end

"""
    nastran_pazy_swept_modal_files(sweep, version; tip=:flexible, nmodes=15) -> Vector{String}

Saves the results of the modal analysis of the swept Technion Pazy wing with `sweep` degrees of sweep (S10, S20), `version` and tip rod `tip` (the .op2 of `nastran_pazy_swept_modal_response`, S<sweep>_<version>_modal<suffix>.op2), to be copied to AeroBeams (test/referenceData/Pazy), in nastran/Pazy_Technion/S<sweep>/<version> and returns their paths: the frequencies of the modes, written by `nastran_pazy_swept_modal_response`, as S<sweep>_<version>_modal_frequencies<suffix>.csv (mode number, eigenvalue [(rad/s)²], frequency [Hz]) and the mode shapes of the first `nmodes` modes as S<sweep>_<version>_modal_shapes<suffix>_wing.csv and S<sweep>_<version>_modal_shapes<suffix>_rod.csv, extracted with extractModeShapes.py (see `nastran_pazy_modal_shapes`; y is the one of the Nastran model, from the root section, and x of the rod is from the leading edge at its y).
"""
function nastran_pazy_swept_modal_files(sweep::Real, version::AbstractString; tip::Symbol=:flexible, nmodes::Integer=15)
    python, nastranDir, caseDir, name = pazy_swept_dirs(sweep, version)
    suffix = pazy_tip_suffix(tip)
    job = "$(name)_modal$suffix"
    shapes = "$(name)_modal_shapes$suffix"

    run(`$python $(joinpath(nastranDir, "extractModeShapes.py")) $caseDir --job $job --nmodes $nmodes --name $shapes`)

    return [joinpath(caseDir, "$(name)_modal_frequencies$suffix.csv"); [joinpath(caseDir, "$(shapes)_$kind.csv") for kind in ("wing", "rod")]]
end

# ------------------------------------------------------------------------------------------------------------
# Technion Pazy swept wings without the tip rod (S10_norod, S20_norod: the decks of pazy_swept_no_rod_bdf)
# ------------------------------------------------------------------------------------------------------------

"""
    nastran_pazy_swept_no_rod_static_response(sweep, version; reference_axis)

Static Nastran workflow for the swept Technion Pazy wing with `sweep` degrees of sweep (S10, S20), `version` = "LE", "TE" or "Baseline", without the tip rod (the deck of `pazy_swept_no_rod_bdf`, which has no property 13, 20, 22 and 24), as `nastran_pazy_swept_static_response` (generatePazy.py --no-rod: the load node is tied to the tip section, with no rod). Files are kept in nastran/Pazy_Technion/S<sweep>/NoRod, named S<sweep>_<version>_norod_RA<x>_static.* (see `pazy_reference_axis_tag`).

Returns `(T, L)`, as in `nastran_pazy_swept_static_response`.

If Nastran is not available on this machine, an error is raised after listing the .bdf to run elsewhere. Copy the .op2 file back to that directory and call this function again. The Python interpreter is the project's `.venv` unless the BEAMROM_PYTHON environment variable is set.
"""
function nastran_pazy_swept_no_rod_static_response(sweep::Real, version::AbstractString; reference_axis::Real)
    python, nastranDir, caseDir, name = pazy_swept_no_rod_dirs(sweep, version)
    bdf = pazy_swept_no_rod_bdf(sweep, version)
    tag = pazy_reference_axis_tag(reference_axis)
    job = "$(name)_$(tag)_static"

    run(`$python $(joinpath(nastranDir, "generatePazy.py")) --bdf $bdf --no-rod --analyses static --sweep $sweep --tag $tag --chord-fraction $reference_axis`)
    run_pazy_nastran(python, nastranDir, caseDir, job)

    L = parse(Float64, match(r"\"L\":\s*([-+0-9.eE]+)", read(joinpath(caseDir, "$job.meta.json"), String)).captures[1])
    T = read_extracted_csv(`$python $(joinpath(nastranDir, "extractResponse.py")) $caseDir --job $job`)

    return T, L
end

"""
    nastran_pazy_swept_no_rod_station_response(sweep, version, y; reference_axis, tol=1e-5, twist_center=:auto)

Generalized displacements of the nodes of a discretized beam of the swept Technion Pazy wing without the tip rod, as `nastran_pazy_swept_station_response` (from the job of `nastran_pazy_swept_no_rod_static_response`, which is generated and run if needed). No file is written.

Returns `(T, L)`, as in `nastran_pazy_swept_station_response`.
"""
function nastran_pazy_swept_no_rod_station_response(sweep::Real, version::AbstractString, y::AbstractVector{<:Real}; reference_axis::Real, tol::Real=1e-5, twist_center::Union{Symbol,Real}=:auto)
    twist_center isa Real || twist_center in (:auto, :axis) || throw(ArgumentError("twist_center must be :auto, :axis or a fraction of the chord, got $twist_center"))
    python, nastranDir, caseDir, name = pazy_swept_no_rod_dirs(sweep, version)
    _, L = nastran_pazy_swept_no_rod_static_response(sweep, version; reference_axis)
    job = "$(name)_$(pazy_reference_axis_tag(reference_axis))_static"

    Tall = read_extracted_csv(`$python $(joinpath(nastranDir, "extractStations.py")) $caseDir --job $job --stations $y --tol $tol --twist-center $twist_center`)
    T = [Tall[6*(k-1)+1:6*k, :] for k in eachindex(y)]

    return T, L
end

"""
    nastran_pazy_swept_no_rod_modal_response(sweep, version; nd=10)

Modal Nastran workflow for the swept Technion Pazy wing without the tip rod (generatePazy.py --no-rod), as `nastran_pazy_swept_modal_response`. Files are kept in nastran/Pazy_Technion/S<sweep>/NoRod, named S<sweep>_<version>_norod_modal*.

Returns `(λ, L)`, as in `nastran_pazy_swept_modal_response`.
"""
function nastran_pazy_swept_no_rod_modal_response(sweep::Real, version::AbstractString; nd::Integer=10)
    python, nastranDir, caseDir, name = pazy_swept_no_rod_dirs(sweep, version)
    bdf = pazy_swept_no_rod_bdf(sweep, version)
    job = "$(name)_modal"
    frequencies = "$(name)_modal_frequencies.csv"

    run(`$python $(joinpath(nastranDir, "generatePazy.py")) --bdf $bdf --no-rod --analyses modal --nd $nd --sweep $sweep`)
    run_pazy_nastran(python, nastranDir, caseDir, job)
    run(`$python $(joinpath(nastranDir, "extractModes.py")) $caseDir --job $job --out $frequencies`)

    L = parse(Float64, match(r"\"L\":\s*([-+0-9.eE]+)", read(joinpath(caseDir, "$job.json"), String)).captures[1])

    return readdlm(joinpath(caseDir, frequencies), ',')[:, 2], L
end

"""
    nastran_pazy_swept_no_rod_modal_files(sweep, version; nmodes=15) -> Vector{String}

Saves the results of the modal analysis of the swept Technion Pazy wing without the tip rod (see `nastran_pazy_swept_no_rod_modal_response`): the frequencies of the modes, S<sweep>_<version>_norod_modal_frequencies.csv, and the mode shapes of the wing of the first `nmodes` modes, S<sweep>_<version>_norod_modal_shapes_wing.csv (extractModeShapes.py writes no rod file, as the model has no rod). Returns their paths.
"""
function nastran_pazy_swept_no_rod_modal_files(sweep::Real, version::AbstractString; nmodes::Integer=15)
    python, nastranDir, caseDir, name = pazy_swept_no_rod_dirs(sweep, version)
    job = "$(name)_modal"
    shapes = "$(name)_modal_shapes"

    run(`$python $(joinpath(nastranDir, "extractModeShapes.py")) $caseDir --job $job --nmodes $nmodes --name $shapes`)

    return [joinpath(caseDir, "$(name)_modal_frequencies.csv"), joinpath(caseDir, "$(shapes)_wing.csv")]
end

# ------------------------------------------------------------------------------------------------------------
# Technion Pazy wing without the tip rod (S0_norod: the same wing and spar, without the tip rod and its tip masses)
# ------------------------------------------------------------------------------------------------------------

const PAZY_NO_ROD_DIR = joinpath("Pazy_Technion", "S0_norod") # folder of the model and of its Nastran files
const PAZY_NO_ROD_NAME = "S0_norod" # name of the base .bdf (bulk data only), prefix of the outputs and of the Nastran jobs

"""
    pazy_no_rod_bdf() -> String

Path of the base .bdf of the Technion Pazy wing without the tip rod, nastran/Pazy_Technion/S0_norod/S0_norod.bdf (it has no rod, its CBEAM property 13 and CONM2 tip masses, nor the RBE2 that ties the rod to the tip).
"""
pazy_no_rod_bdf() = joinpath(@__DIR__, PAZY_NO_ROD_DIR, "$(PAZY_NO_ROD_NAME).bdf")

function pazy_no_rod_dirs()
    python, nastranDir, _ = pazy_python_and_dirs()
    return python, nastranDir, joinpath(nastranDir, PAZY_NO_ROD_DIR)
end

"""
    pazy_no_rod_mass_model_file(; reference_axis, mesh) -> String

Path of the CSV file with the equivalent beam mass model of the Technion Pazy wing without the tip rod for the given reference axis and mesh (see `pazy_beam_mass_model`), which is nastran/Pazy_Technion/S0_norod/S0_norod_RA<x>_inertia_matrices_<mesh>.csv (see `pazy_mesh_tag`).
"""
pazy_no_rod_mass_model_file(; reference_axis::Real, mesh) = joinpath(@__DIR__, PAZY_NO_ROD_DIR, "$(PAZY_NO_ROD_NAME)_$(pazy_reference_axis_tag(reference_axis))_inertia_matrices_$(pazy_mesh_tag(mesh)).csv")

"""
    pazy_no_rod_stiffness_file(; reference_axis) -> String

Path of the CSV file with the global stiffness matrix of the Technion Pazy wing without the tip rod for the given reference axis (chord fraction), which is nastran/Pazy_Technion/S0_norod/S0_norod_RA<x>_stiffness_matrix.csv.
"""
pazy_no_rod_stiffness_file(; reference_axis::Real) = joinpath(@__DIR__, PAZY_NO_ROD_DIR, "$(PAZY_NO_ROD_NAME)_$(pazy_reference_axis_tag(reference_axis))_stiffness_matrix.csv")

"""
    pazy_no_rod_least_squares_stiffness_file(; reference_axis, mesh, segments) -> String

Path of the CSV file with the stiffness matrices of the Technion Pazy wing without the tip rod fitted by `least_squares_stiffness_matrices` (see `pazy_least_squares_stiffness_file`), which is nastran/Pazy_Technion/S0_norod/S0_norod_RA<x>_stiffness_matrices_<mesh>.csv (see `pazy_mesh_tag`), and `_global` after it for `segments = :global` (see `pazy_segments_suffix`).
"""
pazy_no_rod_least_squares_stiffness_file(; reference_axis::Real, mesh, segments::Symbol) = joinpath(@__DIR__, PAZY_NO_ROD_DIR, "$(PAZY_NO_ROD_NAME)_$(pazy_reference_axis_tag(reference_axis))_stiffness_matrices_$(pazy_mesh_tag(mesh))$(pazy_segments_suffix(segments)).csv")

"""
    nastran_pazy_no_rod_static_response(; reference_axis)

Static Nastran workflow for the Technion Pazy wing without the tip rod (nastran/Pazy_Technion/S0_norod/S0_norod.bdf), as `nastran_pazy_static_response` (generatePazy.py --bdf ... --no-rod: the load node is tied to the tip section by an RBE2 of its own, instead of the rod): files are kept in nastran/Pazy_Technion/S0_norod, and each reference axis has its own files, named S0_norod_RA<x>_static.*.

Returns `(T, L)`, as in `nastran_pazy_static_response`.

If Nastran is not available on this machine, an error is raised after listing the .bdf to run elsewhere (e.g. on a Windows VM, with `python nastran/runNastran.py "nastran/Pazy_Technion/S0_norod" --jobs S0_norod_RA0p43_static`). Copy the .op2 file back to that directory and call this function again. The Python interpreter is the project's `.venv` unless the BEAMROM_PYTHON environment variable is set.
"""
function nastran_pazy_no_rod_static_response(; reference_axis::Real)
    python, nastranDir, caseDir = pazy_no_rod_dirs()
    bdf = pazy_no_rod_bdf()
    tag = pazy_reference_axis_tag(reference_axis)
    job = "$(PAZY_NO_ROD_NAME)_$(tag)_static"

    run(`$python $(joinpath(nastranDir, "generatePazy.py")) --bdf $bdf --no-rod --sweep $PAZY_STRAIGHT_SWEEP --analyses static --tag $tag --chord-fraction $reference_axis`)
    run_pazy_nastran(python, nastranDir, caseDir, job)

    L = parse(Float64, match(r"\"L\":\s*([-+0-9.eE]+)", read(joinpath(caseDir, "$job.meta.json"), String)).captures[1])
    T = read_extracted_csv(`$python $(joinpath(nastranDir, "extractResponse.py")) $caseDir --job $job`)

    return T, L
end

"""
    nastran_pazy_no_rod_station_response(y; reference_axis, tol=1e-5, twist_center=:auto)

Generalized displacements of the nodes of a discretized beam of the Technion Pazy wing without the tip rod, under the six unit loads at the end of the reference axis, from the job of `nastran_pazy_no_rod_static_response` (which is generated and run if needed, so no extra Nastran run is required). `y`, `tol` and `twist_center` are as in `nastran_pazy_station_response`, and no file is written.

Returns `(T, L)`, as in `nastran_pazy_station_response`.
"""
function nastran_pazy_no_rod_station_response(y::AbstractVector{<:Real}; reference_axis::Real, tol::Real=1e-5, twist_center::Union{Symbol,Real}=:auto)
    twist_center isa Real || twist_center in (:auto, :axis) || throw(ArgumentError("twist_center must be :auto, :axis or a fraction of the chord, got $twist_center"))
    python, nastranDir, caseDir = pazy_no_rod_dirs()
    _, L = nastran_pazy_no_rod_static_response(; reference_axis)
    job = "$(PAZY_NO_ROD_NAME)_$(pazy_reference_axis_tag(reference_axis))_static"

    Tall = read_extracted_csv(`$python $(joinpath(nastranDir, "extractStations.py")) $caseDir --job $job --stations $y --tol $tol --twist-center $twist_center`)
    T = [Tall[6*(k-1)+1:6*k, :] for k in eachindex(y)]

    return T, L
end

"""
    nastran_pazy_no_rod_modal_response(; nd=10)

Modal Nastran workflow for the Technion Pazy wing without the tip rod (nastran/Pazy_Technion/S0_norod/S0_norod.bdf): generates the modal BDF (generatePazy.py --no-rod, which leaves the tip free, without tie), runs Nastran (runNastran.py) and extracts the eigenvalues (extractModes.py), saved as S0_norod_modal_frequencies.csv. Files are kept in nastran/Pazy_Technion/S0_norod. `nd` is the number of modes requested.

Returns `(λ, L)`, as in `nastran_pazy_modal_response`.

If Nastran is not available on this machine, an error is raised after listing the .bdf to run elsewhere (e.g. on a Windows VM, with `python nastran/runNastran.py "nastran/Pazy_Technion/S0_norod" --jobs S0_norod_modal`). Copy the .op2 file back to that directory and call this function again. The Python interpreter is the project's `.venv` unless the BEAMROM_PYTHON environment variable is set.
"""
function nastran_pazy_no_rod_modal_response(; nd::Integer=10)
    python, nastranDir, caseDir = pazy_no_rod_dirs()
    bdf = pazy_no_rod_bdf()
    job = "$(PAZY_NO_ROD_NAME)_modal"
    frequencies = "$(PAZY_NO_ROD_NAME)_modal_frequencies.csv"

    run(`$python $(joinpath(nastranDir, "generatePazy.py")) --bdf $bdf --no-rod --sweep $PAZY_STRAIGHT_SWEEP --analyses modal --nd $nd`)
    run_pazy_nastran(python, nastranDir, caseDir, job)
    run(`$python $(joinpath(nastranDir, "extractModes.py")) $caseDir --job $job --out $frequencies`)

    L = parse(Float64, match(r"\"L\":\s*([-+0-9.eE]+)", read(joinpath(caseDir, "$job.json"), String)).captures[1])

    return readdlm(joinpath(caseDir, frequencies), ',')[:, 2], L
end

"""
    nastran_pazy_no_rod_modal_files(; nmodes=15) -> Vector{String}

Saves the results of the modal analysis of the Technion Pazy wing without the tip rod (see `nastran_pazy_no_rod_modal_response`): the frequencies of the modes, S0_norod_modal_frequencies.csv, and the mode shapes of the wing of the first `nmodes` modes, S0_norod_modal_shapes_wing.csv (extractModeShapes.py writes no rod file, as the model has no rod). Returns their paths.
"""
function nastran_pazy_no_rod_modal_files(; nmodes::Integer=15)
    python, nastranDir, caseDir = pazy_no_rod_dirs()
    job = "$(PAZY_NO_ROD_NAME)_modal"
    shapes = "$(PAZY_NO_ROD_NAME)_modal_shapes"

    run(`$python $(joinpath(nastranDir, "extractModeShapes.py")) $caseDir --job $job --nmodes $nmodes --name $shapes`)

    return [joinpath(caseDir, "$(PAZY_NO_ROD_NAME)_modal_frequencies.csv"), joinpath(caseDir, "$(shapes)_wing.csv")]
end

# ------------------------------------------------------------------------------------------------------------
# Delft Pazy wing (nastran/Pazy_Delft/pazy_delft.bdf, another FEM of the wing, unswept, like the Technion S0 but without its two ribs close to the root): spar PSHELL 3, tip rod PBEAM 13 tied by RBE2 to the end section of the wing (spar and skin), 17 ribs PSHELL 4 and 10, skin PSHELL 10011 (12 micron, as the Technion S0) and stringers PBEAM 1 and 2 with the roll tied to the spar at the ribs (tieStringerRoll.py --spar-pid 3)
# ------------------------------------------------------------------------------------------------------------

const PAZY_DELFT_DIR = "Pazy_Delft" # folder of the model and of its Nastran files
const PAZY_DELFT_NAME = "pazy_delft" # name of the base .bdf, prefix of the outputs and of the Nastran jobs
const PAZY_DELFT_SPAR_PID = 3
const PAZY_DELFT_ROD_PID = 13
const PAZY_DELFT_RIB_PIDS = (4, 10)
const PAZY_DELFT_SKIP_ROOT = 0 # no ribs close to the root to lump in the root node
const PAZY_DELFT_SKIP_TIP = 3 # the last 3 of the 4 ribs close to the tip are lumped in the tip node (as in the Technion S0)

"""
    pazy_delft_bdf() -> String

Path of the base .bdf (bulk data only) of the Delft Pazy wing, nastran/Pazy_Delft/pazy_delft.bdf.
"""
pazy_delft_bdf() = joinpath(@__DIR__, PAZY_DELFT_DIR, "$(PAZY_DELFT_NAME).bdf")

function pazy_delft_dirs()
    python, nastranDir, _ = pazy_python_and_dirs()
    return python, nastranDir, joinpath(nastranDir, PAZY_DELFT_DIR)
end

"""
    pazy_delft_spar_reference_axis(mesh; reference_axis, sweep=PAZY_STRAIGHT_SWEEP) -> Vector{Vector{Float64}}

Beam nodes along the reference axis of the spar of the Delft Pazy wing (see `spar_reference_axis`), with the spar PSHELL 3 and the ribs PSHELL 4 and 10. The `:nominal` mesh has the root, the first 14 ribs and the tip of the spar (15 elements, the last one with the 4 ribs close to the tip), and `:refined` adds a node in each element but the last.
"""
function pazy_delft_spar_reference_axis(mesh; reference_axis::Real, sweep::Real=PAZY_STRAIGHT_SWEEP)
    return spar_reference_axis(read_gfem(pazy_delft_bdf()), mesh; reference_axis, sweep, spar_pid=PAZY_DELFT_SPAR_PID, rib_pids=PAZY_DELFT_RIB_PIDS, skip_root=PAZY_DELFT_SKIP_ROOT, skip_tip=PAZY_DELFT_SKIP_TIP)
end

"""
    pazy_delft_beam_mass_model(; nodes, reference_axis, sweep=PAZY_STRAIGHT_SWEEP, kwargs...)

Equivalent beam mass model of the Delft Pazy wing (see `pazy_beam_mass_model`) with the nodes of `pazy_delft_spar_reference_axis`.
"""
function pazy_delft_beam_mass_model(; nodes, reference_axis::Real, sweep::Real=PAZY_STRAIGHT_SWEEP, kwargs...)
    return pazy_beam_mass_model(pazy_delft_bdf(); nodes, reference_axis, sweep, spar_pid=PAZY_DELFT_SPAR_PID, rib_pids=PAZY_DELFT_RIB_PIDS, skip_root=PAZY_DELFT_SKIP_ROOT, skip_tip=PAZY_DELFT_SKIP_TIP, kwargs...)
end

"""
    pazy_delft_mass_model_file(; reference_axis, mesh) -> String

Path of the CSV file with the equivalent beam mass model of the Delft Pazy wing for the given reference axis (chord fraction) and mesh, which is nastran/Pazy_Delft/pazy_delft_RA<x>_inertia_matrices_<mesh>.csv (see `pazy_mesh_tag`).
"""
pazy_delft_mass_model_file(; reference_axis::Real, mesh) = joinpath(@__DIR__, PAZY_DELFT_DIR, "$(PAZY_DELFT_NAME)_$(pazy_reference_axis_tag(reference_axis))_inertia_matrices_$(pazy_mesh_tag(mesh)).csv")

"""
    pazy_delft_stiffness_file(; reference_axis, tip=:rigid) -> String

Path of the CSV file with the global stiffness matrix of the Delft Pazy wing for the given reference axis (chord fraction) and tip rod (see `pazy_tip_suffix`), which is nastran/Pazy_Delft/pazy_delft_RA<x>_stiffness_matrix<suffix>.csv.
"""
pazy_delft_stiffness_file(; reference_axis::Real, tip::Symbol=:rigid) = joinpath(@__DIR__, PAZY_DELFT_DIR, "$(PAZY_DELFT_NAME)_$(pazy_reference_axis_tag(reference_axis))_stiffness_matrix$(pazy_tip_suffix(tip)).csv")

"""
    pazy_delft_least_squares_stiffness_file(; reference_axis, mesh, segments, tip=:rigid) -> String

Path of the CSV file with the stiffness matrices of the Delft Pazy wing fitted by `least_squares_stiffness_matrices` (see `pazy_least_squares_stiffness_file`), which is nastran/Pazy_Delft/pazy_delft_RA<x>_stiffness_matrices_<mesh><suffix><segments suffix>.csv.
"""
pazy_delft_least_squares_stiffness_file(; reference_axis::Real, mesh, segments::Symbol, tip::Symbol=:rigid) = joinpath(@__DIR__, PAZY_DELFT_DIR, "$(PAZY_DELFT_NAME)_$(pazy_reference_axis_tag(reference_axis))_stiffness_matrices_$(pazy_mesh_tag(mesh))$(pazy_tip_suffix(tip))$(pazy_segments_suffix(segments)).csv")

"""
    nastran_pazy_delft_static_response(; reference_axis, tip=:rigid)

Static Nastran workflow for the Delft Pazy wing (nastran/Pazy_Delft/pazy_delft.bdf), as `nastran_pazy_static_response` for the Technion wing: generatePazy.py with the spar PSHELL 3 and the tip rod PBEAM 13, runs Nastran (runNastran.py) and extracts the response (extractResponse.py). Files are kept in nastran/Pazy_Delft, named pazy_delft_RA<x>_static<suffix>.* (see `pazy_tip_suffix`).

Returns `(T, L)`, as in `nastran_pazy_static_response`.

If Nastran is not available on this machine, an error is raised after listing the .bdf to run elsewhere (e.g. on a Windows VM, with `python nastran/runNastran.py "nastran/Pazy_Delft" --jobs pazy_delft_RA0p47_static_rigidtip`). Copy the .op2 file back to that directory and call this function again. The Python interpreter is the project's `.venv` unless the BEAMROM_PYTHON environment variable is set.
"""
function nastran_pazy_delft_static_response(; reference_axis::Real, tip::Symbol=:rigid)
    python, nastranDir, caseDir = pazy_delft_dirs()
    tag = pazy_reference_axis_tag(reference_axis)
    job = "$(PAZY_DELFT_NAME)_$(tag)_static" * pazy_tip_suffix(tip)

    run(`$python $(joinpath(nastranDir, "generatePazy.py")) --bdf $(pazy_delft_bdf()) --spar-pid $PAZY_DELFT_SPAR_PID --rod-pid $PAZY_DELFT_ROD_PID --sweep $PAZY_STRAIGHT_SWEEP --analyses static --tag $tag --chord-fraction $reference_axis --tip $tip`)
    run_pazy_nastran(python, nastranDir, caseDir, job)

    L = parse(Float64, match(r"\"L\":\s*([-+0-9.eE]+)", read(joinpath(caseDir, "$job.meta.json"), String)).captures[1])
    T = read_extracted_csv(`$python $(joinpath(nastranDir, "extractResponse.py")) $caseDir --job $job`)

    return T, L
end

"""
    nastran_pazy_delft_station_response(y; reference_axis, tol=1e-5, twist_center=:auto, tip=:rigid)

Generalized displacements of the nodes of a discretized beam of the Delft Pazy wing under the six unit loads at the end of the reference axis, from the job of `nastran_pazy_delft_static_response` (generated and run if needed). `y`, `tol`, `twist_center` and the returned values are as in `nastran_pazy_station_response`; the spar grids at each station are those of PSHELL 3 (extractStations.py --spar-pid), and no file is written.
"""
function nastran_pazy_delft_station_response(y::AbstractVector{<:Real}; reference_axis::Real, tol::Real=1e-5, twist_center::Union{Symbol,Real}=:auto, tip::Symbol=:rigid)
    twist_center isa Real || twist_center in (:auto, :axis) || throw(ArgumentError("twist_center must be :auto, :axis or a fraction of the chord, got $twist_center"))
    python, nastranDir, caseDir = pazy_delft_dirs()
    _, L = nastran_pazy_delft_static_response(; reference_axis, tip)
    job = "$(PAZY_DELFT_NAME)_$(pazy_reference_axis_tag(reference_axis))_static" * pazy_tip_suffix(tip)

    Tall = read_extracted_csv(`$python $(joinpath(nastranDir, "extractStations.py")) $caseDir --job $job --stations $y --tol $tol --twist-center $twist_center --spar-pid $PAZY_DELFT_SPAR_PID`)
    T = [Tall[6*(k-1)+1:6*k, :] for k in eachindex(y)]

    return T, L
end

"""
    nastran_pazy_delft_modal_response(; nd=10, tip=:flexible)

Modal Nastran workflow for the Delft Pazy wing (nastran/Pazy_Delft/pazy_delft.bdf): generates the modal BDF (generatePazy.py, with the tip rod PBEAM 13 rigid or flexible, see `pazy_tip_suffix`), runs Nastran and extracts the eigenvalues (extractModes.py) to pazy_delft_modal_frequencies<suffix>.csv. Files are kept in nastran/Pazy_Delft. `nd` is the number of modes requested.

Returns `(λ, L)`, as in `nastran_pazy_modal_response`.

If Nastran is not available on this machine, an error is raised after listing the .bdf to run elsewhere (e.g. on a Windows VM, with `python nastran/runNastran.py "nastran/Pazy_Delft" --jobs pazy_delft_modal_flextip`). Copy the .op2 file back to that directory and call this function again. The Python interpreter is the project's `.venv` unless the BEAMROM_PYTHON environment variable is set.
"""
function nastran_pazy_delft_modal_response(; nd::Integer=10, tip::Symbol=:flexible)
    python, nastranDir, caseDir = pazy_delft_dirs()
    job = "$(PAZY_DELFT_NAME)_modal" * pazy_tip_suffix(tip)
    frequencies = "$(PAZY_DELFT_NAME)_modal_frequencies$(pazy_tip_suffix(tip)).csv"

    run(`$python $(joinpath(nastranDir, "generatePazy.py")) --bdf $(pazy_delft_bdf()) --spar-pid $PAZY_DELFT_SPAR_PID --rod-pid $PAZY_DELFT_ROD_PID --sweep $PAZY_STRAIGHT_SWEEP --analyses modal --nd $nd --tip $tip`)
    run_pazy_nastran(python, nastranDir, caseDir, job)
    run(`$python $(joinpath(nastranDir, "extractModes.py")) $caseDir --job $job --out $frequencies`)

    L = parse(Float64, match(r"\"L\":\s*([-+0-9.eE]+)", read(joinpath(caseDir, "$job.json"), String)).captures[1])

    return readdlm(joinpath(caseDir, frequencies), ',')[:, 2], L
end

"""
    nastran_pazy_delft_modal_shapes(; tip=:flexible, nmodes=15) -> Vector{String}

Extracts the mode shapes of the first `nmodes` modes of the modal analysis of the Delft Pazy wing with the tip rod `tip` (the .op2 of `nastran_pazy_delft_modal_response`, pazy_delft_modal<suffix>.op2) with extractModeShapes.py, as `nastran_pazy_modal_shapes` does for the Technion wing, with the spar PSHELL 3 and the tip rod PBEAM 13. Saves them as pazy_delft_modal_shapes<suffix>_wing.csv and pazy_delft_modal_shapes<suffix>_rod.csv in nastran/Pazy_Delft (to be copied to AeroBeams, test/referenceData/Pazy) and returns their paths.
"""
function nastran_pazy_delft_modal_shapes(; tip::Symbol=:flexible, nmodes::Integer=15)
    python, nastranDir, caseDir = pazy_delft_dirs()
    suffix = pazy_tip_suffix(tip)
    job = "$(PAZY_DELFT_NAME)_modal" * suffix
    name = "$(PAZY_DELFT_NAME)_modal_shapes" * suffix

    run(`$python $(joinpath(nastranDir, "extractModeShapes.py")) $caseDir --job $job --nmodes $nmodes --name $name --spar-pid $PAZY_DELFT_SPAR_PID --rod-pid $PAZY_DELFT_ROD_PID`)

    return [joinpath(caseDir, "$(name)_$kind.csv") for kind in ("wing", "rod")]
end
