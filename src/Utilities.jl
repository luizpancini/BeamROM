"""
    isotropic_stiffness_matrix(; kwargs...)

Creates a 6x6 sectional stiffness matrix for a cross-section made of isotropic material

# Arguments
- `∞` = value for rigid properties
- `EA` = axial stiffness
- `GAy` = shear stiffness in the x2 direction
- `GAz` = shear stiffness in the x3 direction
- `GJ` = torsional stiffness
- `EIy` = bending stiffness in the x2 direction
- `EIz` = bending stiffness in the x3 direction
- `t2` = offset from reference line to tension center in the x2 direction
- `t3` = offset from reference line to tension center in the x3 direction
- `s2` = offset from reference line to shear center in the x2 direction
- `s3` = offset from reference line to shear center in the x3 direction
"""
function isotropic_stiffness_matrix(; ∞=1e16,EA=∞,GAy=∞,GAz=∞,GJ=∞,EIy=∞,EIz=∞,t2=0,t3=0,s2=0,s3=0)

    # Validate
    @assert all(x->x>0,[EA,GAy,GAz,GJ,EIy,EIz])

    # See Hodges' book Eqs. 4.114 - 4.118
    z = [0 t3 -t2; -s3 0 0; s2 0 0]

    A = diagm([EA,GAy,GAz])

    Tinv = diagm([GJ,EIy,EIz])

    B = A*z

    D = Tinv + z'*A*z

    S = [A B; B' D]

    return S

end
export isotropic_stiffness_matrix