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

# Throws if the stiffness matrix is not symmetric, since only its upper triangle is written
function check_symmetric(S::AbstractMatrix)
    norm(S - S') <= 1e-9 * norm(S) || throw(ArgumentError("the stiffness matrix is not symmetric, so only its upper triangle can be written"))
    return nothing
end

"""
    write_stiffness_csv(path, S)

Writes a symmetric 6x6 sectional stiffness matrix to a CSV file, without header: row i holds the i-th generalized force component (F1, F2, F3, M1, M2, M3) for the six generalized strains from the i-th one on (columns i to 6, the upper triangle, as the lower terms are equal to them), at full precision. The rows have 6, 5, ..., 1 entries.
"""
function write_stiffness_csv(path::AbstractString, S::AbstractMatrix)
    check_symmetric(S)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        for i in axes(S, 1)
            println(io, join((@sprintf("%.14e", S[i, j]) for j in i:size(S, 2)), ","))
        end
    end
    return path
end
export write_stiffness_csv

"""
    write_discretized_stiffness_csv(path, y, S; symmetric=true)

Writes the sectional stiffness matrices `S` of the segments of a discretized beam, with the nodes at the coordinates `y` (`length(S) == length(y) - 1`), to a CSV file with a header: one row per segment with its index, the coordinates of its nodes (`y1`, `y2`) and the entries of its stiffness matrix by row (`S11, S12, ..., S66`, generalized forces F1, F2, F3, M1, M2, M3 by generalized strains), at full precision. With `symmetric=true` (the default) only the upper triangle is written (`S11, S12, ..., S16, S22, ..., S66`), as the matrices must be symmetric; with `symmetric=false` all 36 entries are written, for non-symmetric matrices.
"""
function write_discretized_stiffness_csv(path::AbstractString, y::AbstractVector, S::AbstractVector{<:AbstractMatrix}; symmetric::Bool=true)
    length(S) == length(y) - 1 || throw(ArgumentError("expected $(length(y) - 1) stiffness matrices for $(length(y)) nodes, got $(length(S))"))
    symmetric && foreach(check_symmetric, S)
    entries = [(i, j) for i in 1:6 for j in 1:6 if !symmetric || i <= j]
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "segment,y1,y2,", join(("S$i$j" for (i, j) in entries), ","))
        for k in eachindex(S)
            println(io, k, ",", @sprintf("%.14e", y[k]), ",", @sprintf("%.14e", y[k+1]), ",", join((@sprintf("%.14e", S[k][i, j]) for (i, j) in entries), ","))
        end
    end
    return path
end
export write_discretized_stiffness_csv