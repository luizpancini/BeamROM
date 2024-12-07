"""
    stiffness_matrix(; kwargs...)

Computes the sectional stiffness matrix, given the length of the beam, the loads at each load case and the resulting displacement and rotations

# Keyword arguments
- `F1=1` = force in load case of force in x1-direction
- `F2=1` = force in load case of force in x2-direction
- `F3=1` = force in load case of force in x3-direction
- `M1=1` = moment in load case of moment in x1-direction
- `M2=1` = moment in load case of moment in x2-direction
- `M3=1` = moment in load case of moment in x3-direction
- `L` = length of the beam
- `u_F1` = displacements in load case of force in x1-direction
- `u_F2` = displacements in load case of force in x2-direction
- `u_F3` = displacements in load case of force in x3-direction
- `u_M1` = displacements in load case of moment in x1-direction
- `u_M2` = displacements in load case of moment in x2-direction
- `u_M3` = displacements in load case of moment in x3-direction
- `θ_F1` = rotation in load case of force in x1-direction
- `θ_F2` = rotation in load case of force in x2-direction
- `θ_F3` = rotation in load case of force in x3-direction
- `θ_M1` = rotation in load case of moment in x1-direction
- `θ_M2` = rotation in load case of moment in x2-direction
- `θ_M3` = rotation in load case of moment in x3-direction
"""
function stiffness_matrix(; F1=1,F2=1,F3=1,M1=1,M2=1,M3=1,L,u_F1,u_F2,u_F3,u_M1,u_M2,u_M3,θ_F1,θ_F2,θ_F3,θ_M1,θ_M2,θ_M3)

    # Generalized forces matrix
    F = diagm([F1; F2; F3; M1; M2; M3])
    
    # Generalized displacements matrix
    B = zeros(eltype(u_F1),6,6)
    B[:,1] .= vcat(u_F1,θ_F1)
    B[:,2] .= vcat(u_F2,θ_F2)
    B[:,3] .= vcat(u_F3,θ_F3)
    B[:,4] .= vcat(u_M1,θ_M1)
    B[:,5] .= vcat(u_M2,θ_M2)
    B[:,6] .= vcat(u_M3,θ_M3)

    # Generalized displacement-forces matrix
    D = B * inv(F)

    # Extract submatrices from D
    D11 = D[1:3, 1:3]
    D12 = D[1:3, 4:6]
    D22 = D[4:6, 4:6]

    # Skew-symmetric matrix of the first vector of a triad
    a1tilde = [0 0 0; 0 0 -1; 0 1 0]
    
    # Solve for T, Z, and R using the relations from D components
    T = D22 / L
    Z = (D12 + (L^2/2) * a1tilde * T) / L
    R = (D11 + a1tilde * Z' + (L^3/3) * a1tilde * T * a1tilde - (L^2/2) * Z * a1tilde) / L

    # Construct compliance matrix from R, Z, and T
    C = [R Z; Z' T]
    
    # Invert compliance matrix to find stiffness matrix
    S = inv(C)
    
    return S
end
export stiffness_matrix


# Wrapper function for Jacobian calculation
function stiffness_matrix_wrapper(x,F1,F2,F3,M1,M2,M3,L)

    # Unpack array of generalized displacements
    u_F1 = x[1:3]
    θ_F1 = x[4:6]
    u_F2 = x[7:9]
    θ_F2 = x[10:12]
    u_F3 = x[13:15]
    θ_F3 = x[16:18]
    u_M1 = x[19:21]
    θ_M1 = x[22:24]
    u_M2 = x[25:27]
    θ_M2 = x[28:30]
    u_M3 = x[31:33]
    θ_M3 = x[34:36]

    # Compute stiffness matrix
    S = stiffness_matrix(F1=F1,F2=F2,F3=F3,M1=M1,M2=M2,M3=M3,L=L,u_F1=u_F1,u_F2=u_F2,u_F3=u_F3,u_M1=u_M1,u_M2=u_M2,u_M3=u_M3,θ_F1=θ_F1,θ_F2=θ_F2,θ_F3=θ_F3,θ_M1=θ_M1,θ_M2=θ_M2,θ_M3=θ_M3)

    return vcat(S...)
end


"""
    stiffness_matrix_derivatives(; kwargs...)

Computes the derivatives (Jacobian terms) of each of the entries of the sectional stiffness matrix with respect to the generalized displacements resulting from the loads

# Keyword arguments
- `F1=1` = force in load case of force in x1-direction
- `F2=1` = force in load case of force in x2-direction
- `F3=1` = force in load case of force in x3-direction
- `M1=1` = moment in load case of moment in x1-direction
- `M2=1` = moment in load case of moment in x2-direction
- `M3=1` = moment in load case of moment in x3-direction
- `L` = length of the beam
- `u_F1` = displacements in load case of force in x1-direction
- `u_F2` = displacements in load case of force in x2-direction
- `u_F3` = displacements in load case of force in x3-direction
- `u_M1` = displacements in load case of moment in x1-direction
- `u_M2` = displacements in load case of moment in x2-direction
- `u_M3` = displacements in load case of moment in x3-direction
- `θ_F1` = rotation in load case of force in x1-direction
- `θ_F2` = rotation in load case of force in x2-direction
- `θ_F3` = rotation in load case of force in x3-direction
- `θ_M1` = rotation in load case of moment in x1-direction
- `θ_M2` = rotation in load case of moment in x2-direction
- `θ_M3` = rotation in load case of moment in x3-direction
"""
function stiffness_matrix_derivatives(; F1=1,F2=1,F3=1,M1=1,M2=1,M3=1,L,u_F1,u_F2,u_F3,u_M1,u_M2,u_M3,θ_F1,θ_F2,θ_F3,θ_M1,θ_M2,θ_M3)

    # Array of generalized displacements
    genDisp = vcat(u_F1,θ_F1,u_F2,θ_F2,u_F3,θ_F3,u_M1,θ_M1,u_M2,θ_M2,u_M3,θ_M3)

    # Jacobian matrix
    J = ForwardDiff.jacobian(x -> stiffness_matrix_wrapper(x,F1,F2,F3,M1,M2,M3,L), genDisp)

    # Jacobian terms
    # S_11
    S11_uF1,S11_θF1,S11_uF2,S11_θF2,S11_uF3,S11_θF3,S11_uM1,S11_θM1,S11_uM2,S11_θM2,S11_uM3,S11_θM3 = J[1,1:3],J[1,4:6],J[1,7:9],J[1,10:12],J[1,13:15],J[1,16:18],J[1,19:21],J[1,22:24],J[1,25:27],J[1,28:30],J[1,31:33],J[1,34:36]
    # S_12
    S12_uF1,S12_θF1,S12_uF2,S12_θF2,S12_uF3,S12_θF3,S12_uM1,S12_θM1,S12_uM2,S12_θM2,S12_uM3,S12_θM3 = J[2,1:3],J[2,4:6],J[2,7:9],J[2,10:12],J[2,13:15],J[2,16:18],J[2,19:21],J[2,22:24],J[2,25:27],J[2,28:30],J[2,31:33],J[2,34:36]
    # S_13
    S13_uF1,S13_θF1,S13_uF2,S13_θF2,S13_uF3,S13_θF3,S13_uM1,S13_θM1,S13_uM2,S13_θM2,S13_uM3,S13_θM3 = J[3,1:3],J[3,4:6],J[3,7:9],J[3,10:12],J[3,13:15],J[3,16:18],J[3,19:21],J[3,22:24],J[3,25:27],J[3,28:30],J[3,31:33],J[3,34:36]
    # S_14
    S14_uF1,S14_θF1,S14_uF2,S14_θF2,S14_uF3,S14_θF3,S14_uM1,S14_θM1,S14_uM2,S14_θM2,S14_uM3,S14_θM3 = J[4,1:3],J[4,4:6],J[4,7:9],J[4,10:12],J[4,13:15],J[4,16:18],J[4,19:21],J[4,22:24],J[4,25:27],J[4,28:30],J[4,31:33],J[4,34:36]
    # S_15
    S15_uF1,S15_θF1,S15_uF2,S15_θF2,S15_uF3,S15_θF3,S15_uM1,S15_θM1,S15_uM2,S15_θM2,S15_uM3,S15_θM3 = J[5,1:3],J[5,4:6],J[5,7:9],J[5,10:12],J[5,13:15],J[5,16:18],J[5,19:21],J[5,22:24],J[5,25:27],J[5,28:30],J[5,31:33],J[5,34:36]
    # S_16
    S16_uF1,S16_θF1,S16_uF2,S16_θF2,S16_uF3,S16_θF3,S16_uM1,S16_θM1,S16_uM2,S16_θM2,S16_uM3,S16_θM3 = J[6,1:3],J[6,4:6],J[6,7:9],J[6,10:12],J[6,13:15],J[6,16:18],J[6,19:21],J[6,22:24],J[6,25:27],J[6,28:30],J[6,31:33],J[6,34:36]
    # S_22
    S22_uF1,S22_θF1,S22_uF2,S22_θF2,S22_uF3,S22_θF3,S22_uM1,S22_θM1,S22_uM2,S22_θM2,S22_uM3,S22_θM3 = J[8,1:3],J[8,4:6],J[8,7:9],J[8,10:12],J[8,13:15],J[8,16:18],J[8,19:21],J[8,22:24],J[8,25:27],J[8,28:30],J[8,31:33],J[8,34:36]
    # S_23
    S23_uF1,S23_θF1,S23_uF2,S23_θF2,S23_uF3,S23_θF3,S23_uM1,S23_θM1,S23_uM2,S23_θM2,S23_uM3,S23_θM3 = J[9,1:3],J[9,4:6],J[9,7:9],J[9,10:12],J[9,13:15],J[9,16:18],J[9,19:21],J[9,22:24],J[9,25:27],J[9,28:30],J[9,31:33],J[9,34:36]
    # S_24
    S24_uF1,S24_θF1,S24_uF2,S24_θF2,S24_uF3,S24_θF3,S24_uM1,S24_θM1,S24_uM2,S24_θM2,S24_uM3,S24_θM3 = J[10,1:3],J[10,4:6],J[10,7:9],J[10,10:12],J[10,13:15],J[10,16:18],J[10,19:21],J[10,22:24],J[10,25:27],J[10,28:30],J[10,31:33],J[10,34:36]
    # S_25
    S25_uF1,S25_θF1,S25_uF2,S25_θF2,S25_uF3,S25_θF3,S25_uM1,S25_θM1,S25_uM2,S25_θM2,S25_uM3,S25_θM3 = J[11,1:3],J[11,4:6],J[11,7:9],J[11,10:12],J[11,13:15],J[11,16:18],J[11,19:21],J[11,22:24],J[11,25:27],J[11,28:30],J[11,31:33],J[11,34:36]
    # S_26
    S26_uF1,S26_θF1,S26_uF2,S26_θF2,S26_uF3,S26_θF3,S26_uM1,S26_θM1,S26_uM2,S26_θM2,S26_uM3,S26_θM3 = J[12,1:3],J[12,4:6],J[12,7:9],J[12,10:12],J[12,13:15],J[12,16:18],J[12,19:21],J[12,22:24],J[12,25:27],J[12,28:30],J[12,31:33],J[12,34:36]
    # S_33
    S33_uF1,S33_θF1,S33_uF2,S33_θF2,S33_uF3,S33_θF3,S33_uM1,S33_θM1,S33_uM2,S33_θM2,S33_uM3,S33_θM3 = J[15,1:3],J[15,4:6],J[15,7:9],J[15,10:12],J[15,13:15],J[15,16:18],J[15,19:21],J[15,22:24],J[15,25:27],J[15,28:30],J[15,31:33],J[15,34:36]
    # S_34
    S34_uF1,S34_θF1,S34_uF2,S34_θF2,S34_uF3,S34_θF3,S34_uM1,S34_θM1,S34_uM2,S34_θM2,S34_uM3,S34_θM3 = J[16,1:3],J[16,4:6],J[16,7:9],J[16,10:12],J[16,13:15],J[16,16:18],J[16,19:21],J[16,22:24],J[16,25:27],J[16,28:30],J[16,31:33],J[16,34:36]
    # S_35
    S35_uF1,S35_θF1,S35_uF2,S35_θF2,S35_uF3,S35_θF3,S35_uM1,S35_θM1,S35_uM2,S35_θM2,S35_uM3,S35_θM3 = J[17,1:3],J[17,4:6],J[17,7:9],J[17,10:12],J[17,13:15],J[17,16:18],J[17,19:21],J[17,22:24],J[17,25:27],J[17,28:30],J[17,31:33],J[17,34:36]
    # S_36
    S36_uF1,S36_θF1,S36_uF2,S36_θF2,S36_uF3,S36_θF3,S36_uM1,S36_θM1,S36_uM2,S36_θM2,S36_uM3,S36_θM3 = J[18,1:3],J[18,4:6],J[18,7:9],J[18,10:12],J[18,13:15],J[18,16:18],J[18,19:21],J[18,22:24],J[18,25:27],J[18,28:30],J[18,31:33],J[18,34:36]
    # S_44
    S44_uF1,S44_θF1,S44_uF2,S44_θF2,S44_uF3,S44_θF3,S44_uM1,S44_θM1,S44_uM2,S44_θM2,S44_uM3,S44_θM3 = J[22,1:3],J[22,4:6],J[22,7:9],J[22,10:12],J[22,13:15],J[22,16:18],J[22,19:21],J[22,22:24],J[22,25:27],J[22,28:30],J[22,31:33],J[22,34:36]
    # S_45
    S45_uF1,S45_θF1,S45_uF2,S45_θF2,S45_uF3,S45_θF3,S45_uM1,S45_θM1,S45_uM2,S45_θM2,S45_uM3,S45_θM3 = J[23,1:3],J[23,4:6],J[23,7:9],J[23,10:12],J[23,13:15],J[23,16:18],J[23,19:21],J[23,22:24],J[23,25:27],J[23,28:30],J[23,31:33],J[23,34:36]
    # S_46
    S46_uF1,S46_θF1,S46_uF2,S46_θF2,S46_uF3,S46_θF3,S46_uM1,S46_θM1,S46_uM2,S46_θM2,S46_uM3,S46_θM3 = J[24,1:3],J[24,4:6],J[24,7:9],J[24,10:12],J[24,13:15],J[24,16:18],J[24,19:21],J[24,22:24],J[24,25:27],J[24,28:30],J[24,31:33],J[24,34:36]
    # S_55
    S55_uF1,S55_θF1,S55_uF2,S55_θF2,S55_uF3,S55_θF3,S55_uM1,S55_θM1,S55_uM2,S55_θM2,S55_uM3,S55_θM3 = J[29,1:3],J[29,4:6],J[29,7:9],J[29,10:12],J[29,13:15],J[29,16:18],J[29,19:21],J[29,22:24],J[29,25:27],J[29,28:30],J[29,31:33],J[29,34:36]
    # S_56
    S56_uF1,S56_θF1,S56_uF2,S56_θF2,S56_uF3,S56_θF3,S56_uM1,S56_θM1,S56_uM2,S56_θM2,S56_uM3,S56_θM3 = J[30,1:3],J[30,4:6],J[30,7:9],J[30,10:12],J[30,13:15],J[30,16:18],J[30,19:21],J[30,22:24],J[30,25:27],J[30,28:30],J[30,31:33],J[30,34:36]
    # S_66
    S66_uF1,S66_θF1,S66_uF2,S66_θF2,S66_uF3,S66_θF3,S66_uM1,S66_θM1,S66_uM2,S66_θM2,S66_uM3,S66_θM3 = J[36,1:3],J[36,4:6],J[36,7:9],J[36,10:12],J[36,13:15],J[36,16:18],J[36,19:21],J[36,22:24],J[36,25:27],J[36,28:30],J[36,31:33],J[36,34:36]

    return J,S11_uF1,S11_θF1,S11_uF2,S11_θF2,S11_uF3,S11_θF3,S11_uM1,S11_θM1,S11_uM2,S11_θM2,S11_uM3,S11_θM3,S12_uF1,S12_θF1,S12_uF2,S12_θF2,S12_uF3,S12_θF3,S12_uM1,S12_θM1,S12_uM2,S12_θM2,S12_uM3,S12_θM3,S13_uF1,S13_θF1,S13_uF2,S13_θF2,S13_uF3,S13_θF3,S13_uM1,S13_θM1,S13_uM2,S13_θM2,S13_uM3,S13_θM3,S14_uF1,S14_θF1,S14_uF2,S14_θF2,S14_uF3,S14_θF3,S14_uM1,S14_θM1,S14_uM2,S14_θM2,S14_uM3,S14_θM3,S15_uF1,S15_θF1,S15_uF2,S15_θF2,S15_uF3,S15_θF3,S15_uM1,S15_θM1,S15_uM2,S15_θM2,S15_uM3,S15_θM3,S16_uF1,S16_θF1,S16_uF2,S16_θF2,S16_uF3,S16_θF3,S16_uM1,S16_θM1,S16_uM2,S16_θM2,S16_uM3,S16_θM3,S22_uF1,S22_θF1,S22_uF2,S22_θF2,S22_uF3,S22_θF3,S22_uM1,S22_θM1,S22_uM2,S22_θM2,S22_uM3,S22_θM3,S23_uF1,S23_θF1,S23_uF2,S23_θF2,S23_uF3,S23_θF3,S23_uM1,S23_θM1,S23_uM2,S23_θM2,S23_uM3,S23_θM3,S24_uF1,S24_θF1,S24_uF2,S24_θF2,S24_uF3,S24_θF3,S24_uM1,S24_θM1,S24_uM2,S24_θM2,S24_uM3,S24_θM3,S25_uF1,S25_θF1,S25_uF2,S25_θF2,S25_uF3,S25_θF3,S25_uM1,S25_θM1,S25_uM2,S25_θM2,S25_uM3,S25_θM3,S26_uF1,S26_θF1,S26_uF2,S26_θF2,S26_uF3,S26_θF3,S26_uM1,S26_θM1,S26_uM2,S26_θM2,S26_uM3,S26_θM3,S33_uF1,S33_θF1,S33_uF2,S33_θF2,S33_uF3,S33_θF3,S33_uM1,S33_θM1,S33_uM2,S33_θM2,S33_uM3,S33_θM3,S34_uF1,S34_θF1,S34_uF2,S34_θF2,S34_uF3,S34_θF3,S34_uM1,S34_θM1,S34_uM2,S34_θM2,S34_uM3,S34_θM3,S35_uF1,S35_θF1,S35_uF2,S35_θF2,S35_uF3,S35_θF3,S35_uM1,S35_θM1,S35_uM2,S35_θM2,S35_uM3,S35_θM3,S36_uF1,S36_θF1,S36_uF2,S36_θF2,S36_uF3,S36_θF3,S36_uM1,S36_θM1,S36_uM2,S36_θM2,S36_uM3,S36_θM3,S44_uF1,S44_θF1,S44_uF2,S44_θF2,S44_uF3,S44_θF3,S44_uM1,S44_θM1,S44_uM2,S44_θM2,S44_uM3,S44_θM3,S45_uF1,S45_θF1,S45_uF2,S45_θF2,S45_uF3,S45_θF3,S45_uM1,S45_θM1,S45_uM2,S45_θM2,S45_uM3,S45_θM3,S46_uF1,S46_θF1,S46_uF2,S46_θF2,S46_uF3,S46_θF3,S46_uM1,S46_θM1,S46_uM2,S46_θM2,S46_uM3,S46_θM3,S55_uF1,S55_θF1,S55_uF2,S55_θF2,S55_uF3,S55_θF3,S55_uM1,S55_θM1,S55_uM2,S55_θM2,S55_uM3,S55_θM3,S56_uF1,S56_θF1,S56_uF2,S56_θF2,S56_uF3,S56_θF3,S56_uM1,S56_θM1,S56_uM2,S56_θM2,S56_uM3,S56_θM3,S66_uF1,S66_θF1,S66_uF2,S66_θF2,S66_uF3,S66_θF3,S66_uM1,S66_θM1,S66_uM2,S66_θM2,S66_uM3,S66_θM3
end
export stiffness_matrix_derivatives


"""
    displacements_from_loads(; kwargs...)

Computes the generalized displacements for a cantilevered, prismatic beam with given length and sectional stiffness matrix, under some specified load case

# Keyword arguments
- `loadCase` = load case (unit force or moment in some direction)
- `L` = length of the beam
- `S` = sectional stiffness matrix
"""
function displacements_from_loads(; loadCase, L, S)

    # Set generalized forces array
    if loadCase == "all"
        F = Matrix(1.0*LinearAlgebra.I,6,6)
    else
        F1=F2=F3=M1=M2=M3 = 0
        if loadCase == "F1"
            F1 = 1
        elseif loadCase == "F2"
            F2 = 1
        elseif loadCase == "F3"
            F3 = 1
        elseif loadCase == "M1"   
            M1 = 1
        elseif loadCase == "M2"   
            M2 = 1
        elseif loadCase == "M3"   
            M3 = 1
        end
        F = [F1; F2; F3; M1; M2; M3]
    end

    # Compliance matrix and submatrices
    C = inv(S)
    R = C[1:3,1:3]
    Z = C[1:3,4:6]
    T = C[4:6,4:6]

    # Skew-symmetric matrix of the first vector of a triad
    a1tilde = [0 0 0; 0 0 -1; 0 1 0]

    # Generalized displacement-force matrix
    D11 = L*R + L^2/2*Z*a1tilde-a1tilde*Z' - L^3/3*a1tilde*T*a1tilde
    D12 = L*Z - L^2/2*a1tilde*T
    D21 = L*Z' + L^2/2*T*a1tilde
    D22 = L*T

    D = [D11 D12; D21 D22]

    # Generalized displacements
    uθ = D * F

    return uθ
end
export displacements_from_loads