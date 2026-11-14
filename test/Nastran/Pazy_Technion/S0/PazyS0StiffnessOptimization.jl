using BeamROM, AeroBeams, DelimitedFiles, LinearAlgebra, Metaheuristics, Printf, Random
include(joinpath(pkgdir(BeamROM), "nastran", "workflow.jl"))

# Heuristic optimization (Metaheuristics.jl) of the stiffness matrix(es) of the Pazy wing equivalent beam, so that the first natural frequencies of AeroBeams match those of the Nastran modal analysis (SOL 103) of the full FEM.
# The model is the same as in PazyS0AeroBeamsModal.jl (clamped beam, no gravity), with the least squares stiffness matrices of PazyS0LeastSquaresStiffnessModel.jl (those of the chosen segments) as the nominal matrices (with GAy and GAz modified), which are the reference and the starting point of the search.
# Run PazyS0MassModel.jl, PazyS0LeastSquaresStiffnessModel.jl (with the same segments) and PazyS0ModalAnalysis.jl first if the CSV files are missing or outdated. Run this one in an environment with AeroBeams and Metaheuristics.

# Optimization method -- choose one of :ABC (Artificial Bee Colony), :DE (Differential Evolution), :ECA (Evolutionary Centers Algorithm), :GA (Genetic Algorithm), :PSO (Particle Swarm Optimization), :SA (Simulated Annealing)
method = :PSO

# Reference axis and beam mesh of the mass and stiffness models to read (as in PazyS0MassModel.jl), which select the files of their combination. The optimized results are saved in files that also depend on them
# mesh must be :nominal or :refined
reference_axis = 0.4653
mesh = :nominal

# Tip rod of the Nastran model: :rigid (the rod is rigid, so the tip section of the wing is, as the stiffness matrices require) or :flexible (the rod of the model, which bends). The results of each have the suffix _rigidtip or _flextip in their file names (see pazy_tip_suffix)
tip = :rigid

# Segments (as in PazyS0LeastSquaresStiffnessModel.jl), each with its own stiffness matrix and design variables, shared by its elements:
# :nodes = one segment for each element (between consecutive beam nodes)
# :global = a single segment for the whole beam
# :ends = one segment for each of the two elements of the end segments (root and tip) and a single one for the middle segment, for a total of five segments
segments = :nodes

# Starting point of the search: :nominal = the least squares stiffness matrices of the segments (S0_RA<RA>_stiffness_matrices_<mesh>mesh_<tip>tip[_global|_ends].csv), or :best = the best solution so far with the same segments (S0_RA<RA>_best_stiffness_matrices_<mesh>mesh[_global|_ends].csv, which must exist). With :best, the matrix of each segment is the mean of the element matrices of that segment in the best solution, and the bounds are centered on it. Note that the best file does not record GAy, GAz and the nominal matrices: use the same ones as in the run that generated it
startPoint = :nominal

# Shear stiffnesses [N] that replace those of the nominal matrices of all the segments (useful for segments=:global). Use nothing to keep the original values (they must be positive)
GAy = nothing
GAz = nothing

# Number of modes in the modal objective function
modalWeight = 1
nModes = 5

# Flutter criterion (in addition to the modal one): relative errors of the flutter onset and offset speeds of mode `flutterMode` with respect to the reference speeds, for each root pitch angle. The objective is: (sum of the squared relative frequency errors of the first nModes modes) + flutterSpeedWeight * (mean of the squared relative errors of the available reference speeds). Use flutterSpeedWeight = 0 to disable the criterion (the aeroelastic analyses are skipped if flutterStrainsWeight = 0 too)
flutterSpeedWeight = 1
flutterMode = 3
flutterRootAngles = [3, 5, 7]            # [deg]
flutterOnsetRef = [49, 43, 38]          # [m/s]
flutterOffsetRef = [55, 48, 43]    # [m/s] (missing = not considered)
aeroSolver = Indicial()
airfoil = FlatPlate
# Airspeeds range for flutter computation
flutterURange = vcat(36:2:64)
# Number of tracked oscillatory modes (must be at least flutterMode) 
nModesFlutter = 3
# Squared relative error assigned to a reference speed (or strain) if the mode does not flutter (onset) within flutterURange
flutterMissingError = 1.0

# 1/6-span strains criterion (in addition to the modal and flutter ones): relative errors of the axial strains [microns] at the quarter-span of the wing, at the flutter onset and offset of mode `flutterMode`, for the root pitch angles of flutterRootAngles. The strains are those at the 95 mm spanwise position strain gauge: ϵ11 - κ2 * x3SG, with x3SG on the top of the spar. The objective is then: (modal error) + flutterSpeedWeight * (flutter error) + flutterStrainsWeight * (mean of the squared relative errors of the available reference strains). Use flutterStrainsWeight = 0 to disable the criterion
flutterStrainsWeight = 0.25
strainsOnsetRef = [missing, 1350, 1450]     # [microns] (missing = not considered)
strainsOffsetRef = [missing, 1900, 1850]    # [microns] (missing = not considered)
x2SG, x3SG = 50e-3, 2.25e-3 / 2   # [m] (top of LE of the spar, whose cross-section is 60 x 2.25 mm)

# Design variables of each matrix (indices of the generalized forces/strains: 1 = F1, 2 = F2, 3 = F3, 4 = M1, 5 = M2, 6 = M3).
# Diagonal terms: S[i,i] = S0[i,i] * exp(x), with x in [x_start - log(diagonalRange[i]), x_start + log(diagonalRange[i])] (diagonalRange has one factor ≥ 1 for each of the 6 diagonal terms), which keeps them positive and symmetric by construction
# Off-diagonal terms: S[i,j] = S[j,i] = S0[i,j] + x * sqrt(S0[i,i] * S0[j,j]) for the free ones (freeCouplings), with x ∈ [x_start - couplingRange[k], x_start + couplingRange[k]] (couplingRange has one value for each of the free couplings), so S[i,j] varies by up to couplingRange[k] * sqrt(S0[i,i] * S0[j,j]) around its starting value. The others keep the nominal correlation coefficient ρ0 with the current diagonal terms
# The bounds are relative to the starting point (nominal or best), so they are centered on it: an out-of-bounds best solution is not an error. Matrices that are not positive definite (possible within the bounds if couplingRange is large) are penalized in the objective function
freeDiagonals = [4]
diagonalRange = [1, 1, 1, 1.5, 1, 1]
freeCouplings = [(4,5)]
couplingRange = [1e-1]

# Options for the optimization algorithm
@assert method in [:ABC, :DE, :ECA, :GA, :PSO, :SA] "selected method must be one of [:ABC, :DE, :ECA, :GA, :PSO, :SA]"
if method == :ABC
    seed = 1
    colonySize = 50
    limit = 10
    maxIterations = 100
elseif method == :DE
    seed = 1
    populationSize = 50
    maxIterations = 100
    F = 0.5
    CR = 0.9
elseif method == :ECA
    seed = 1
    N = 0
    maxIterations = 100
elseif method == :GA
    seed = 1
    populationSize = 50
    p_mutation = 1e-4
    p_crossover = 0.5
    maxIterations = 100
elseif method == :PSO
    seed = rand(Int)
    seed < 0 && (seed *= -1)
    swarmSize = 50
    C1 = 2.0
    C2 = 2.0
    ω = 0.8
    maxIterations = 100
elseif method == :SA
    seed = 1
    N = 500
    tol_fun = 1e-4
    maxIterations = 20
end

# Files
mesh in (:nominal, :refined) || error("mesh must be :nominal or :refined")
segments in (:nodes, :global, :ends) || error("segments must be :nodes, :global or :ends")
dir = normpath(pkgdir(BeamROM), "nastran", "Pazy_Technion", "S0")
massFile = pazy_mass_model_file(; reference_axis, mesh)
stiffnessFile = pazy_least_squares_stiffness_file(; reference_axis, mesh, segments, tip)
modalFile = joinpath(dir, "S0_modal_frequencies$(pazy_tip_suffix(tip)).csv")
caseTag = "$(pazy_reference_axis_tag(reference_axis))"
caseSuffix = "$(pazy_mesh_tag(mesh))$(pazy_segments_suffix(segments))"
outputFile = joinpath(dir, "S0_$(caseTag)_optimized_stiffness_matrices_$caseSuffix.csv")

# Inertia data (as in PazyS0AeroBeamsModal.jl)
mass = readdlm(massFile, ',', Float64; skipstart=1)
nNodes = size(mass, 1)
nElem = nNodes - 1
y = mass[:, 3]
L = y[end] - y[1]
m = mass[:, 5]
offsets = [mass[k, 6:8] for k in 1:nNodes]
tensors = [[mass[k, 9]  mass[k, 12] mass[k, 13];
            mass[k, 12] mass[k, 10] mass[k, 14];
            mass[k, 13] mass[k, 14] mass[k, 11]] for k in 1:nNodes]
R = [0.0 1.0 0.0; 
    -1.0 0.0 0.0; 
     0.0 0.0 1.0]
offsets = [R * o for o in offsets]
tensors = [R * J * R' for J in tensors]

# Sectional inertia matrices, with each interior node split in halves between its two elements
tilde(η) = [0 -η[3] η[2]; η[3] 0 -η[1]; -η[2] η[1] 0]
sectional_inertia(m, η, J, Δℓ) = 1 / Δℓ * [m * I(3) -m * tilde(η); m * tilde(η) J - m * tilde(η)^2]
I_sec = [zeros(6, 6) for _ in 1:nElem]
for k in 1:nElem
    Δℓ = y[k+1] - y[k]
    wStart = k == 1 ? 1.0 : 0.5
    wEnd = k == nElem ? 1.0 : 0.5
    I_sec[k] += sectional_inertia(wStart * m[k], [-Δℓ / 2, 0, 0] + offsets[k], wStart * tensors[k], Δℓ)
    I_sec[k] += sectional_inertia(wEnd * m[k+1], [Δℓ / 2, 0, 0] + offsets[k+1], wEnd * tensors[k+1], Δℓ)
end

# Segments: the elements of each segment share the same stiffness matrix (the same as in PazyS0LeastSquaresStiffnessModel.jl)
if segments == :nodes
    elementSegment = collect(1:nElem)
elseif segments == :ends
    elementSegment = [1; 2; fill(3, nElem - 4); 4; 5]
else
    elementSegment = ones(Int, nElem)
end
nSeg = maximum(elementSegment)
breakNodes = [findfirst(==(g), elementSegment) for g in 1:nSeg]
breaks = [y[breakNodes] .- y[1]; y[end] - y[1]]

# Symmetric 6x6 matrix from its upper triangle by rows (S11, S12, ..., S16, S22, ..., S66), as saved by write_discretized_stiffness_csv
function symmetric_from_upper(v)
    S = zeros(6, 6)
    k = 0
    for i in 1:6, j in i:6
        k += 1
        S[i, j] = S[j, i] = v[k]
    end
    return S
end

# Nominal stiffness matrices of the segments (least squares), symmetrized, with modified shear stiffnesses
lsData = readdlm(stiffnessFile, ',', Float64; skipstart=1)
size(lsData, 1) == nSeg || error("$stiffnessFile has $(size(lsData, 1)) segments, but the mass model has $nSeg: run PazyS0LeastSquaresStiffnessModel.jl with the same reference_axis, mesh and segments")
isapprox(lsData[:, 2:3], [breaks[1:end-1] breaks[2:end]]; atol=1e-6 * L) || error("the segments of $stiffnessFile do not match those of the mass model: run PazyS0LeastSquaresStiffnessModel.jl with the same reference_axis, mesh and segments")
S0s = [symmetric_from_upper(lsData[g, 4:end]) for g in 1:nSeg]
for S0 in S0s
    S0 .= (S0 + S0') / 2
    for (i, GA) in ((2, GAy), (3, GAz))
        isnothing(GA) && continue
        S0[i, :] .= 0
        S0[:, i] .= 0
        S0[i, i] = GA
    end
    all(diag(S0) .> 0) || error("the diagonal terms of the nominal stiffness matrices must be positive: set GAy and GAz")
end
modalWeight >= 0 || error("modalWeight must be non-negative")
flutterSpeedWeight >= 0 || error("flutterSpeedWeight must be non-negative")
flutterStrainsWeight >= 0 || error("flutterStrainsWeight must be non-negative")
length(flutterOnsetRef) == length(flutterOffsetRef) == length(flutterRootAngles) || error("flutterRootAngles, flutterOnsetRef and flutterOffsetRef must have the same length")
length(strainsOnsetRef) == length(strainsOffsetRef) == length(flutterRootAngles) || error("flutterRootAngles, strainsOnsetRef and strainsOffsetRef must have the same length")
1 <= flutterMode <= nModesFlutter || error("flutterMode must be between 1 and nModesFlutter")
ρ0s = [S0 ./ (sqrt.(diag(S0)) * sqrt.(diag(S0))') for S0 in S0s]

# Nastran frequencies [Hz], from the eigenvalues [(rad/s)²]
λNastran = readdlm(modalFile, ',', Float64)[:, 2]
length(λNastran) >= nModes || error("the Nastran modal analysis has only $(length(λNastran)) modes: increase nd in PazyS0ModalAnalysis.jl")
fNastran = sqrt.(λNastran[1:nModes]) / (2π)

# Mapping from the design variables x to the stiffness matrices. The variables are divided in one block per segment: in each block, first the free diagonal terms (log of the scale factors), then the free couplings (normalized increments of the off-diagonal terms, in units of sqrt(S0[i,i] * S0[j,j]), with S0 the nominal matrix of the segment)
nD = length(freeDiagonals)
nBlock = nD + length(freeCouplings)
nVars = nSeg * nBlock
function block_matrix(x, g)
    S0, ρ0 = S0s[g], ρ0s[g]
    d = copy(diag(S0))
    for (k, i) in enumerate(freeDiagonals)
        d[i] *= exp(x[k])
    end
    S = ρ0 .* (sqrt.(d) * sqrt.(d)')
    for (k, (i, j)) in enumerate(freeCouplings)
        S[i, j] = S[j, i] = S0[i, j] + x[nD+k] * sqrt(S0[i, i] * S0[j, j])
    end
    ρ = S ./ (sqrt.(d) * sqrt.(d)')
    return Symmetric(S), ρ
end
function stiffness_matrices(x)
    blocks = [block_matrix(x[(g-1)*nBlock+1:g*nBlock], g) for g in 1:nSeg]
    return [blocks[elementSegment[k]][1] for k in 1:nElem], [b[1] for b in blocks], [b[2] for b in blocks]
end

# Natural frequencies [Hz] of AeroBeams for the stiffness matrices of the elements
function frequencies(Selem)
    beam = create_Beam(name="wing", length=L, nElements=nElem, normalizedNodalPositions=(y .- y[1]) ./ L, S=[Matrix(Sk) for Sk in Selem], I=I_sec)
    clamp = create_BC(name="clamp", beam=beam, node=1, conditions=Dict(:u1A=>0, :u2A=>0, :u3A=>0, :p1A=>0, :p2A=>0, :p3A=>0))
    model = create_Model(name="Pazy", beams=[beam], BCs=[clamp], gravityVector=[0, 0, 0])
    problem = create_EigenProblem(model=model, eigenmethod=:krylovkit, nev=40)
    solve!(problem)
    return sort([imag(λ) / (2π) for λ in problem.eigensolution.λ if imag(λ) > 1e-6])
end

# Objective function: sum of the squared relative frequency errors, with a penalty if any matrix is not positive definite or if the eigenproblem cannot be solved
penalty = 1e3
function frequency_error(Selem)
    f = try
        frequencies(Selem)
    catch
        return 2 * penalty
    end
    length(f) >= nModes || return 2 * penalty
    return sum(((f[1:nModes] .- fNastran) ./ fNastran) .^ 2)
end

# Flutter onset and offset speeds [m/s] of mode flutterMode and root strains [microns] at them (missing if not found within flutterURange) for the stiffness matrices of the elements and the root pitch angle θ [rad]. The aeroelastic model is the same as in create_Pazy (AeroBeams) except for the stiffness and inertia matrices and the nodes, which are those of the optimization
_, _, chord, normSparPos, _ = geometrical_properties_Pazy()
function flutter_results(Selem, θ)
    surf = create_AeroSurface(ppSolver=aeroSolver, airfoil=deepcopy(airfoil), c=chord, normSparPos=normSparPos, hasTipCorrection=true, tipLossFunction=U -> tip_loss_function_Pazy("Exponential", θ, U, 0.0), tipLossFunctionIsAirspeedDependent=true)
    beam = create_Beam(name="wingSpar", length=L, nElements=nElem, normalizedNodalPositions=(y .- y[1]) ./ L, S=[Matrix(Sk) for Sk in Selem], I=I_sec, p0=transform_to_WM([-π / 2, 0, θ]; from=:E231), surface=surf)
    clamp = create_BC(name="clamp", beam=beam, node=1, conditions=Dict(:u1A=>0, :u2A=>0, :u3A=>0, :p1A=>0, :p2A=>0, :p3A=>0))
    model = create_Model(name="Pazy", beams=[beam], BCs=[clamp], gravityVector=[0, 0, -standard_atmosphere(0).g])
    nU = length(flutterURange)
    problems = Vector{EigenProblem}(undef, nU)
    eigensols = Vector{EigenSolution}(undef, nU)
    oneSixthSpanEps = zeros(nU)
    for (i, U) in enumerate(flutterURange)
        update_v_A!(model, v_A=[0, U, 0])
        problems[i] = create_EigenProblem(model=model, eigenmethod=:krylovkit, nev=40, x0=i == 1 ? zeros(0) : problems[i-1].systemArrays.x)
        solve!(problems[i])
        eigensols[i] = problems[i].eigensolution
        oneSixthSpanState = problems[i].loadstepHistory.states[end][3].midpoint
        oneSixthSpanEps[i] = oneSixthSpanState.γ[1] - oneSixthSpanState.κ[2] * x3SG - oneSixthSpanState.κ[3] * x2SG
    end
    freqs, damps, modeShapes, _ = mode_tracking(eigensols)
    rads2Hz!(freqs)
    freqsOsc, dampsOsc, _ = filter_tracked_modes(freqs, damps, modeShapes; fmin=0.1, maxTracks=nModesFlutter)
    onsetSpeeds, _, onsetModes, onsetEps, offsetSpeeds, _, offsetModes, offsetEps = find_flutter(flutterURange, freqsOsc, dampsOsc, oneSixthSpanEps)
    iOnset = findfirst(==(flutterMode), onsetModes)
    isnothing(iOnset) && return missing, missing, missing, missing
    onset = onsetSpeeds[iOnset]
    offsets = offsetSpeeds[offsetModes .== flutterMode]
    offsetStrains = offsetEps[offsetModes .== flutterMode]
    iOffset = findfirst(>(onset), offsets)
    return onset, isnothing(iOffset) ? missing : offsets[iOffset], 1e6 * onsetEps[iOnset], isnothing(iOffset) ? missing : 1e6 * offsetStrains[iOffset]
end

# Mean of the squared relative errors of the flutter onset and offset speeds, and mean of those of the root strains at them, with respect to the reference ones, for the root pitch angles with references. If the mode does not flutter within the range, the error is flutterMissingError. If it has no offset within the range, the speed is taken as the end of the range, and the strain error is flutterMissingError. The mean is 0 if there are no references. The aeroelastic analysis is skipped for the angles without references of the criteria with non-zero weight
function flutter_errors(Selem)
    flutterErrors = Float64[]
    strainsErrors = Float64[]
    rel_error(v, vref) = ismissing(v) ? flutterMissingError : ((v - vref) / vref)^2
    for (k, θ) in enumerate(flutterRootAngles)
        useFlutter = flutterSpeedWeight > 0 && !(ismissing(flutterOnsetRef[k]) && ismissing(flutterOffsetRef[k]))
        useStrains = flutterStrainsWeight > 0 && !(ismissing(strainsOnsetRef[k]) && ismissing(strainsOffsetRef[k]))
        (useFlutter || useStrains) || continue
        onset, offset, onsetStrain, offsetStrain = flutter_results(Selem, deg2rad(θ))
        if useFlutter
            ismissing(flutterOnsetRef[k]) || push!(flutterErrors, rel_error(onset, flutterOnsetRef[k]))
            ismissing(flutterOffsetRef[k]) || push!(flutterErrors, rel_error(ismissing(onset) ? missing : coalesce(offset, flutterURange[end]), flutterOffsetRef[k]))
        end
        if useStrains
            ismissing(strainsOnsetRef[k]) || push!(strainsErrors, rel_error(onsetStrain, strainsOnsetRef[k]))
            ismissing(strainsOffsetRef[k]) || push!(strainsErrors, rel_error(offsetStrain, strainsOffsetRef[k]))
        end
    end
    mean_or_zero(v) = isempty(v) ? 0.0 : sum(v) / length(v)
    return mean_or_zero(flutterErrors), mean_or_zero(strainsErrors)
end

# Total error: modal error + flutterSpeedWeight * flutter error + flutterStrainsWeight * strains error (the aeroelastic analyses are skipped if the modal problem fails or if both weights are 0)
function total_error(Selem)
    modal = frequency_error(Selem)
    (flutterSpeedWeight == 0 && flutterStrainsWeight == 0) && return modalWeight * modal
    modal >= penalty && return modalWeight * modal
    flutter, strains = try
        flutter_errors(Selem)
    catch
        return 2 * penalty
    end
    return modalWeight * modal + flutterSpeedWeight * flutter + flutterStrainsWeight * strains
end
function objective(x)
    Selem, _, ρs = stiffness_matrices(x)
    for ρ in ρs
        any(abs.(ρ[i, j] for i in 1:6, j in 1:6 if i != j) .>= 1) && return penalty
        λmin = eigmin(Symmetric(ρ))
        λmin > 1e-6 || return penalty * (1 + 1e-6 - λmin)
    end
    return total_error(Selem)
end

# Table with one row per element (segment, y1, y2 and the 36 terms of the stiffness matrix by rows), in the format of the discretized stiffness matrices (PazyS0DiscretizedStiffnessModel.jl)
function element_table(Selem)
    header = ["segment", "y1", "y2", ["S$i$j" for i in 1:6 for j in 1:6]...]
    rows = [[k, y[k] - y[1], y[k+1] - y[1], vec(Matrix(Selem[k])')...] for k in 1:nElem]
    return [permutedims(header); permutedims.(rows)...]
end

# Best solution found so far with these segments, shared by all runs with the same reference axis, mesh and segments: the stiffness matrices of all the elements and the value of the objective function, which depends only on these matrices, on nModes, on modalWeight and on the flutter and strains criteria settings. The file is updated only by a run that finds a lower objective. Its first two lines are "objective,nModes,flutterKey,modalWeight" (flutterKey identifies the flutter and strains settings, and is 0 if flutterSpeedWeight = 0 and flutterStrainsWeight = 0) and the values, followed by the element table
bestFile = joinpath(dir, "S0_$(caseTag)_best_stiffness_matrices_$caseSuffix.csv")
flutterKey = if flutterSpeedWeight == 0 && flutterStrainsWeight == 0
    UInt64(0)
else
    flutterSettings = (flutterSpeedWeight, flutterMode, flutterRootAngles, flutterOnsetRef, flutterOffsetRef, flutterURange, nModesFlutter, flutterMissingError, string(aeroSolver), airfoil.name)
    flutterStrainsWeight == 0 ? hash(flutterSettings) : hash((flutterSettings, flutterStrainsWeight, strainsOnsetRef, strainsOffsetRef))
end
function write_best(file, Selem, obj)
    open(file, "w") do io
        println(io, "objective,nModes,flutterKey,modalWeight")
        println(io, repr(obj), ",", nModes, ",", flutterKey, ",", modalWeight)
        writedlm(io, element_table(Selem), ',')
    end
end
function read_best(file)
    isfile(file) || return nothing
    values = split(readlines(file)[2], ',')
    obj = parse(Float64, values[1])
    nModesBest = parse(Int, values[2])
    flutterKeyBest = length(values) >= 3 ? parse(UInt64, values[3]) : UInt64(0)
    # Files without modalWeight were written with the unit weight
    modalWeightBest = length(values) >= 4 ? parse(Float64, values[4]) : 1.0
    data = readdlm(file, ',', Float64; skipstart=3)
    size(data, 1) == nElem || error("$file has $(size(data, 1)) elements, but the mass model has $nElem: remove it if the mesh was changed")
    Selem = [Symmetric(Matrix(transpose(reshape(data[k, 4:end], 6, 6)))) for k in 1:nElem]
    # The stored objective is not comparable if it was computed with another number of modes, modal weight or flutter or strains settings
    (nModesBest == nModes && flutterKeyBest == flutterKey && modalWeightBest == modalWeight) || (obj = total_error(Selem))
    return (objective=obj, S=Selem)
end
best = read_best(bestFile)

# Bounds
length(diagonalRange) == 6 && all(diagonalRange .>= 1) || error("diagonalRange must have 6 entries, all greater than or equal to 1")
couplingRange = collect(Float64, couplingRange)
length(couplingRange) == length(freeCouplings) && all(couplingRange .>= 0) || error("couplingRange must have one non-negative value for each of the free couplings")
widthBlock = [log.(diagonalRange[freeDiagonals]); couplingRange]

# Design variables of the stiffness matrices of the elements (inverse of block_matrix): for each segment, the log scale factors of the free diagonal terms and the normalized increments of the free couplings, from the mean of the element matrices of the segment
function design_variables(Selem)
    xs = zeros(nVars)
    for g in 1:nSeg
        S0 = S0s[g]
        diagS0 = diag(S0)
        elems = findall(==(g), elementSegment)
        Sg = sum(Matrix(Selem[k]) for k in elems) / length(elems)
        dg = diag(Sg)
        sg = sqrt.(dg)
        o = (g - 1) * nBlock
        for (k, i) in enumerate(freeDiagonals)
            xs[o+k] = log(dg[i] / diagS0[i])
        end
        for (k, (i, j)) in enumerate(freeCouplings)
            xs[o+nD+k] = (Sg[i, j] - S0[i, j]) / sqrt(diagS0[i] * diagS0[j])
        end
    end
    return xs
end

# Nominal solution
xNominal = zeros(nVars)
S_nominal = [Symmetric(S0) for S0 in S0s]
f_nominal = frequencies(stiffness_matrices(xNominal)[1])[1:nModes]
obj0 = objective(xNominal)
@printf("Nominal objective: %.4e\n", obj0)
isnothing(best) || @printf("Best objective so far: %.4e (%s)\n", best.objective, basename(bestFile))

# Starting point of the search
startPoint in (:nominal, :best) || error("startPoint must be :nominal or :best")
if startPoint == :best
    isnothing(best) && error("startPoint = :best requires $bestFile, which does not exist: run the optimization with startPoint = :nominal first")
    x0 = design_variables(best.S)
    @printf("Starting from the best solution so far: objective %.4e\n", objective(x0))
else
    x0 = xNominal
end

# Bounds of the design variables, relative to the starting point
lower = x0 .- repeat(widthBlock, nSeg)
upper = x0 .+ repeat(widthBlock, nSeg)
bounds = [lower'; upper']

# Search
Random.seed!(seed)
algorithm = if method == :ABC
    ABC(N=colonySize, limit=limit, options=Options(iterations=maxIterations, seed=seed, parallel_evaluation=false, verbose=true))
elseif method == :DE 
    DE(N=populationSize, F=F, CR=CR, options=Options(iterations=maxIterations, seed=seed, parallel_evaluation=false, verbose=true))
elseif method == :ECA
    ECA(N=N, options=Options(iterations=maxIterations, seed=seed, parallel_evaluation=false, verbose=true))
elseif method == :GA
    GA(N=populationSize, p_crossover=p_crossover, p_mutation=p_mutation, mutation=SlightMutation(), options=Options(iterations=maxIterations, seed=seed, parallel_evaluation=false, verbose=true))
elseif method == :PSO
    PSO(N=swarmSize, C1=C1, C2=C2, ω=ω, options=Options(iterations=maxIterations, seed=seed, parallel_evaluation=false, verbose=true))
elseif method == :SA
    SA(N=N, tol_fun=tol_fun, options=Options(iterations=maxIterations, seed=seed, parallel_evaluation=false, verbose=true))
end
set_user_solutions!(algorithm, x0, objective)
result = optimize(objective, bounds, algorithm)
x = minimizer(result)
S_elem, S_seg, ρ_seg = stiffness_matrices(x)
f_opt = frequencies(S_elem)[1:nModes]

# Checks of the constraints
for g in 1:nSeg
    issymmetric(S_seg[g]) || error("the optimized matrix of segment $g is not symmetric")
    all(diag(S_seg[g]) .> 0) || error("the optimized matrix of segment $g has non-positive diagonal terms")
end

# Results
@printf("\nStiffness matrices: %s (%d segments, %d design variables), starting from %s\n", segments, nSeg, length(x), startPoint)
@printf("Objective: %.4e (nominal) -> %.4e (optimized), %d function evaluations\n", obj0, minimum(result), result.f_calls)
@printf("%-4s %14s %14s %14s %12s %12s\n", "Mode", "Nastran [Hz]", "Nominal [Hz]", "Optimized [Hz]", "Nominal [%]", "Optimized [%]")
for i in 1:nModes
    @printf("%-4d %14.4f %14.4f %14.4f %12.2f %12.2f\n", i, fNastran[i], f_nominal[i], f_opt[i], 100 * (f_nominal[i] - fNastran[i]) / fNastran[i], 100 * (f_opt[i] - fNastran[i]) / fNastran[i])
end
if flutterSpeedWeight > 0 || flutterStrainsWeight > 0
    fmt(v) = ismissing(v) ? @sprintf("%6s", "-") : @sprintf("%6.1f", v)
    fmtEps(v) = ismissing(v) ? @sprintf("%6s", "-") : @sprintf("%6.0f", v)
    S_nominal_elem = stiffness_matrices(xNominal)[1]
    nominalResults = [flutter_results(S_nominal_elem, deg2rad(θ)) for θ in flutterRootAngles]
    optimizedResults = [flutter_results(S_elem, deg2rad(θ)) for θ in flutterRootAngles]
    @printf("\nFlutter speeds of mode %d [m/s] (reference / nominal / optimized)\n", flutterMode)
    @printf("%-8s %-24s %-24s\n", "θ [deg]", "Onset", "Offset")
    for (k, θ) in enumerate(flutterRootAngles)
        nom, opt = nominalResults[k], optimizedResults[k]
        @printf("%-8g %s / %s / %s    %s / %s / %s\n", θ, fmt(flutterOnsetRef[k]), fmt(nom[1]), fmt(opt[1]), fmt(flutterOffsetRef[k]), fmt(nom[2]), fmt(opt[2]))
    end
    @printf("\nRoot strains at the flutter speeds of mode %d [microns] (reference / nominal / optimized)\n", flutterMode)
    @printf("%-8s %-24s %-24s\n", "θ [deg]", "Onset", "Offset")
    for (k, θ) in enumerate(flutterRootAngles)
        nom, opt = nominalResults[k], optimizedResults[k]
        @printf("%-8g %s / %s / %s    %s / %s / %s\n", θ, fmtEps(strainsOnsetRef[k]), fmtEps(nom[3]), fmtEps(opt[3]), fmtEps(strainsOffsetRef[k]), fmtEps(nom[4]), fmtEps(opt[4]))
    end
end
for g in 1:nSeg
    elems = findall(==(g), elementSegment)
    println(nSeg == 1 ? "\nDiagonal terms of the stiffness matrix (nominal -> optimized):" : "\nSegment $g (elements $(join(elems, ", "))), diagonal terms of the stiffness matrix (nominal -> optimized):")
    for i in 1:6
        @printf("  S[%d,%d] = %12.4e -> %12.4e (x%.3f)\n", i, i, S_nominal[g][i, i], S_seg[g][i, i], S_seg[g][i, i] / S_nominal[g][i, i])
    end
end
println("\nOptimized stiffness matri" * (nSeg == 1 ? "x:" : "ces:"))
for g in 1:nSeg
    nSeg > 1 && println("Segment $g:")
    display(Matrix(S_seg[g]))
end

# Write the optimized matrices of the segments, in the same format as the least squares ones (PazyS0LeastSquaresStiffnessModel.jl)
write_discretized_stiffness_csv(outputFile, breaks, [Matrix(Sg) for Sg in S_seg])
println("\nOptimized stiffness matrices written to $outputFile")

# Update the best solution, if improved
objOpt = minimum(result)
if isnothing(best) || objOpt < best.objective
    write_best(bestFile, S_elem, objOpt)
    @printf("New best solution: objective %.4e (previous: %s), written to %s\n", objOpt, isnothing(best) ? "none" : @sprintf("%.4e", best.objective), bestFile)
else
    @printf("The best solution so far (objective %.4e) was not improved by this run (%.4e)\n", best.objective, objOpt)
end
