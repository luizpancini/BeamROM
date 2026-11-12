using BeamROM, Plots, LinearAlgebra, ForwardDiff

# ==============================================================================
# User input data
# ==============================================================================
# Generalized displacement and stiffness entry to relate
displacement = :u_F3       # Choose :u_F1, :θ_M2, etc.
stiffnessEntry = :S_33     # Choose :S_11, :S_22, etc.

# Variation range for displacement (±Δvar) and number of points
Δvar = 2e-3
nPoints = 400+1

# Beam length
L = 1.0

# Generalized displacements due to unit loads
disp_data = Dict(
    :u_F1 => [4.996872249e-8; 0; 0],
    :u_F2 => [0; 8.008386067e-5; 0],
    :u_F3 => [0; 0; 4.973927885e-2],
    :u_M1 => [0; 0; -3.003244423e-8],
    :u_M2 => [0; 0; -7.474333048e-2],
    :u_M3 => [0; 1.19941913e-4; 0],
    :θ_F1 => [0.0, 0.0, 0.0],
    :θ_F2 => [0; 0; 1.19931791e-4],
    :θ_F3 => [0; -7.473818123e-2; 0],
    :θ_M1 => [9.787431948e-2; 7.881497155e-7; 0],
    :θ_M2 => [0; 1.4955157e-1; 0],
    :θ_M3 => [0; 0; 2.39458926e-4],
)

# ==============================================================================
# Validation and indices setting
# ==============================================================================
# Indices of stiffness entry
idx_i = parse(Int, string(stiffnessEntry)[3])
idx_j = parse(Int, string(stiffnessEntry)[4])
if !(idx_i in [1,2,3]) || !(idx_j in [1,2,3])
    error("Set stiffnessEntry as 'S_ij', i and j between 1 and 3")
end

# Non-zero index of chosen displacement vector
indicesMap = Dict(
    :u_F1 => 1, :u_F2 => 2, :u_F3 => 3,
    :u_M1 => 3, :u_M2 => 3, :u_M3 => 2,
    :θ_F1 => 1, :θ_F2 => 3, :θ_F3 => 2,
    :θ_M1 => 1, :θ_M2 => 2, :θ_M3 => 3
)
idx_vec = indicesMap[displacement]

# Jacobian column and row for corresponding displacement
dispOrder = [:u_F1, :θ_F1, :u_F2, :θ_F2, :u_F3, :θ_F3, 
             :u_M1, :θ_M1, :u_M2, :θ_M2, :u_M3, :θ_M3]
pos = findfirst(==(displacement), dispOrder)
if pos === nothing
    error("Invalid displacement")
end
idx_jacobian_col = (pos - 1) * 3 + idx_vec
idx_jacobian_row = idx_i + 6 * (idx_j - 1)

# ==============================================================================
# Expected stiffness matrix (S0) and corresponding chosen entry (S0_nominal)
# ==============================================================================
b, h = 50e-3, 2e-3
E, G = 200e9, 76.92e9
A = b * h
Iy, Iz = (b * h^3) / 12, (h * b^3) / 12
J_torsion = Iy + Iz
Kt = 6.5e-3
Ksy = Ksz = 5/6

S0 = diagm([E*A, G*A*Ksy, G*A*Ksz, G*J_torsion*Kt, E*Iy, E*Iz])
S0_nominal = S0[idx_i, idx_j]

# ==============================================================================
# Analysis
# ==============================================================================
# Nominal and range of values for analysis
nominalValue = disp_data[displacement][idx_vec]
valueRange = range(nominalValue * (1 - Δvar), stop=nominalValue * (1 + Δvar), length=nPoints)

# Initialize arrays
S_values = Float64[]
J_abs_values = Float64[]
J_norm_values = Float64[]

# Loop displacement values
for v in valueRange
    # Current displacement value
    disp_temp = deepcopy(disp_data)
    disp_temp[displacement][idx_vec] = v
    
    # Keyword arguments for calculations
    kwargs = [k => v for (k, v) in disp_temp]
    
    # Stiffness matrix with current displacement value
    S = stiffness_matrix(; L=L, kwargs...)

    # Absolute and normalized Jacobians
    J_norm, J_abs = stiffness_matrix_normalized_derivatives(; L=L, kwargs...)
    
    # Set arrays
    push!(S_values, S[idx_i, idx_j])
    push!(J_abs_values, J_abs[idx_jacobian_row, idx_jacobian_col])
    push!(J_norm_values, J_norm[idx_jacobian_row, idx_jacobian_col])
end

# ==============================================================================
# Plotting
# ==============================================================================
# Settings
label_x = "$(displacement)[$(idx_vec)]"
lw = 2

# Stiffness entry x displacement value
p1 = plot(valueRange, S_values,
          ylims=[max(minimum(S_values),-5*S0_nominal),Inf],
          label="Computed", 
          linewidth=lw, 
          color=:blue,
          xlabel=label_x, 
          ylabel="$(stiffnessEntry)", 
          title="Numerical sensitivity analysis: $(stiffnessEntry) vs $(label_x)")
hline!(p1, [S0_nominal], label="Expected", linestyle=:dash, color=:red, linewidth=lw)

# Normalized Jacobian (∂S/∂u * u/S)
p2 = plot(valueRange, J_norm_values, 
          label="Normalized Jacobian", 
          color=:green, 
          linewidth=lw,
          xlabel=label_x,
          ylabel="∂$(stiffnessEntry)/∂$(label_x) $(label_x)/$(stiffnessEntry) ")

# Absolute Jacobian (∂S/∂u)
p3 = plot(valueRange, J_abs_values, 
          label="Absolute Jacobian", 
          color=:orange, 
          linewidth=lw,
          xlabel=label_x,
          ylabel="∂$(stiffnessEntry)/∂$(label_x)")

# Subplots
plot(p1, p2, p3, layout=(3, 1), size=(900, 900), margin=5Plots.mm)