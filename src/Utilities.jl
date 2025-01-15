"""
    isotropic_stiffness_matrix(; kwargs...)

Creates a 6x6 sectional stiffness matrix for a cross-section made of isotropic material

# Arguments
- `∞` = value for rigid properties
- `EA` = axial stiffness
- `GA2` = shear stiffness in the x2 direction
- `GA3` = shear stiffness in the x3 direction
- `GJ` = torsional stiffness
- `EI2` = bending stiffness in the x2 direction
- `EI3` = bending stiffness in the x3 direction
- `t2` = offset from reference line to tension center in the x2 direction
- `t3` = offset from reference line to tension center in the x3 direction
- `s2` = offset from reference line to shear center in the x2 direction
- `s3` = offset from reference line to shear center in the x3 direction
"""
function isotropic_stiffness_matrix(; ∞=1e16,EA=∞,GA2=∞,GA3=∞,GJ=∞,EI2=∞,EI3=∞,t2=0,t3=0,s2=0,s3=0)

    # Validate
    @assert all(x->x>0,[EA,GA2,GA3,GJ,EI2,EI3])

    # See Hodges' book Eqs. 4.114 - 4.118
    z = [0 t3 -t2; -s3 0 0; s2 0 0]

    A = diagm([EA,GA2,GA3])

    Tinv = diagm([GJ,EI2,EI3])

    B = A*z

    D = Tinv + z'*A*z

    S = [A B; B' D]

    return S
end
export isotropic_stiffness_matrix


"""
    anisotropic_stiffness_matrix(; kwargs...)

Creates a 6x6 sectional stiffness matrix for a cross-section made of anisotropic material

# Arguments
- `∞` = value for rigid properties
- `EA` = axial stiffness
- `GA2` = shear stiffness in the x2 direction
- `GA3` = shear stiffness in the x3 direction
- `GJ` = torsional stiffness
- `EI2` = bending stiffness in the x2 direction
- `EI3` = bending stiffness in the x3 direction
- `S12` = axial-shear (in x2 direction) coupling
- `S13` = axial-shear (in x3 direction) coupling
- `S14` = axial-torsion coupling
- `S15` = axial-bending (in x2 direction) coupling
- `S16` = axial-bending (in x3 direction) coupling
- `S23` = shear-shear coupling
- `S24` = shear (in x2 direction) - torsion coupling
- `S25` = shear-bending (in x2 direction) coupling
- `S26` = shear (in x2 direction) - bending (in x3 direction) coupling
- `S34` = shear (in x3 direction) - torsion coupling
- `S35` = shear (in x3 direction) - bending (in x2 direction) coupling
- `S36` = shear-bending (in x3 direction) coupling
- `S45` = torsion-bending (in x2 direction) coupling
- `S46` = torsion-bending (in x3 direction) coupling
- `S56` = bending-bending coupling
"""
function anisotropic_stiffness_matrix(; ∞=1e16,EA=∞,GA2=∞,GA3=∞,GJ=∞,EI2=∞,EI3=∞,S12=0,S13=0,S14=0,S15=0,S16=0,S23=0,S24=0,S25=0,S26=0,S34=0,S35=0,S36=0,S45=0,S46=0,S56=0)

    # Validate
    @assert all(x->x>0,[EA,GA2,GA3,GJ,EI2,EI3])

    S = [EA S12 S13 S14 S15 S16;
        S12 GA2 S23 S24 S25 S26;
        S13 S23 GA3 S34 S35 S36;
        S14 S24 S34  GJ S45 S46;
        S15 S25 S35 S45 EI2 S56;
        S16 S26 S36 S46 S56 EI3]

    return S
end
export anisotropic_stiffness_matrix