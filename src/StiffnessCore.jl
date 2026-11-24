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
    D = B/F

    return stiffness_matrix_from_D(D, L)
end

"""
    stiffness_matrix_from_D(D, L)

Computes the sectional stiffness matrix from the generalized displacement-force matrix `D` (6 x 6, rows [u; θ] and columns [F1, F2, F3, M1, M2, M3]) at the free end of a prismatic cantilever beam of length `L`

For a prismatic elastic beam `D` is symmetric (reciprocity: `D21 = D12'`), but only `D11`, `D12` and `D22` are used, so with data that are not exactly reciprocal the resulting compliance, and hence the stiffness matrix, are not symmetric (`D21` is ignored, and `R` and `T` are taken as they are)
"""
function stiffness_matrix_from_D(D, L)

    # Extract submatrices from D
    D11 = D[1:3, 1:3]
    D12 = D[1:3, 4:6]
    D22 = D[4:6, 4:6]

    # Skew-symmetric matrix of the first vector of a triad
    a1tilde = [0 0 0; 0 0 -1; 0 1 0]
    
    # Solve for T, Z, and R using the relations from D components
    T = D22 / L
    Z = (D12 + (L^2/2) * a1tilde * T) / L
    R = (D11 + (L^2/2) * a1tilde * Z' - (L^2/2) * Z * a1tilde + (L^3/3) * a1tilde * T * a1tilde) / L

    # Construct compliance matrix from R, Z, and T
    C = [R Z; Z' T]
        
    # Invert compliance matrix to find stiffness matrix
    S = inv(C)
    
    return S
end

"""
    discretized_stiffness_matrices(; y, T, kwargs...)

Computes the sectional stiffness matrices of the `n-1` segments between the `n` nodes of a beam (discretized stiffness condensation), clamped at the first node, from the generalized displacements of the nodes under the six load cases applied at the last node

Each segment is treated as a prismatic cantilever clamped at its upstream node and loaded at its downstream node by the generalized forces transmitted through it, which are the applied force and the applied moment plus the moment of that force about the downstream node. The generalized displacements of the segment are those of the downstream node relative to the upstream one, in the frame that follows the rigid body rotation of the upstream node (`remove_rigid_rotation=true`), and then `stiffness_matrix` applies to the segment. Returns a vector with the (6 x 6) stiffness matrix of each segment, from the root

# Keyword arguments
- `y` = coordinates of the nodes along the beam axis (x1-direction), increasing from the clamped node
- `T` = vector with the (6 x 6) generalized displacements of each node, rows [u; θ] and columns [F1, F2, F3, M1, M2, M3]; the clamped node has zeros
- `y_load=y[end]` = coordinate of the point where the loads are applied
- `F1=1`, `F2=1`, `F3=1`, `M1=1`, `M2=1`, `M3=1` = loads in each load case
"""
function discretized_stiffness_matrices(; y, T, y_load=y[end], F1=1, F2=1, F3=1, M1=1, M2=1, M3=1)

    n = length(y)
    @assert n >= 2 "at least two nodes are required"
    @assert length(T) == n "T must have the generalized displacements of each node"
    @assert all(diff(y) .> 0) "y must be increasing"
    @assert all(size(Tk) == (6, 6) for Tk in T) "the generalized displacements of each node must be a 6 x 6 matrix"

    # Skew-symmetric matrix of the first vector of a triad
    a1tilde = [0 0 0; 0 0 -1; 0 1 0]
    I3 = Matrix(1.0*LinearAlgebra.I, 3, 3)

    # Generalized forces matrix
    F = diagm([F1; F2; F3; M1; M2; M3])

    # Initialize array to store stiffness matrices of each segment
    S = Vector{Matrix{Float64}}(undef, n - 1)

    # Loop over each segment
    for k in 1:n-1
        # Length of the current segment
        Δℓ = y[k+1] - y[k]

        # Generalized displacements of the downstream node relative to the upstream one: subtract from the relative displacement of a segment the displacement of its downstream node due to the rotation of its upstream node (the segment is clamped at its upstream end)
        ΔB = Matrix{Float64}(T[k+1]) - Matrix{Float64}(T[k])
        ΔB[1:3, :] .+= Δℓ * a1tilde * T[k][4:6, :]

        # Generalized forces at the downstream node: force F and moment M + (y_load - y_{k+1}) e1 × F
        Φ = [I3 zeros(3, 3); (y_load - y[k+1]) * a1tilde I3] * F

        # Stiffness matrix of the current segment
        S[k] = stiffness_matrix_from_D(ΔB/Φ, Δℓ)
    end

    return S
end

"""
    clamped_segment_stiffness_matrices(; y, T, intervals=1, method=:average, weights=nothing, kwargs...)

Computes the sectional stiffness matrices of the `n-1` segments between the `n` nodes of a beam from the analysis of windows of `intervals` consecutive segments as cantilevers: the section of the upstream node of a window is clamped and the six load cases are applied at its downstream node, where the generalized displacements are taken. Since the loads are applied where the displacements are taken, the displacement-force matrix of each window is reciprocal (symmetric) if the model is, unlike that of `discretized_stiffness_matrices`, which takes the relative displacements of a single analysis (loads at the end of the beam). Returns a vector with the (6 x 6) stiffness matrix of each segment, from the root

The stiffness matrix of a window is that of a prismatic cantilever of its length. With `intervals > 1` the end effects of the clamp and of the load introduction (warping), which stiffen a short beam, are spread over a longer one. The windows start at every node (`n - intervals` windows, window `w` having the segments `w` to `w + intervals - 1`), and the stiffness matrix of a segment is a combination of those of the windows that contain it (see `method`), so the segments at the ends of the beam have fewer windows. For `intervals = 1` each segment is its own window

# Keyword arguments
- `y` = coordinates of the nodes along the beam axis (x1-direction), increasing from the root
- `T` = vector with the `n - intervals` (6 x 6) generalized displacements of the downstream node of each window, rows [u; θ] and columns [F1, F2, F3, M1, M2, M3]
- `intervals=1` = number of segments of each window
- `method=:average` = how the stiffness matrices of the windows that contain a segment are combined: `:average` is their mean and `:weighted_average` is their mean weighted by the position of the segment in each window (`weights`), normalized by the sum of the weights of the windows that contain it
- `weights=nothing` = for `:weighted_average`, vector of `intervals` weights of the segment by its position in the window, from the one next to the clamped node to the one next to the loaded node. Default `1:intervals`, that is, proportional to the distance to the clamp: the clamp restrains the warping and stiffens the window, mostly next to it, so the windows in which the segment is farther from the clamp weigh more. With the Technion Pazy wing, this brought the torsional stiffness S44 of the segments closer to the value of its uniform part (windows of 4 intervals) than the mean, and the reverse weights (`intervals:-1:1`) were farther. It is not allowed with `:average`
- `F1=1`, `F2=1`, `F3=1`, `M1=1`, `M2=1`, `M3=1` = loads in each load case
"""
function clamped_segment_stiffness_matrices(; y, T, intervals=1, method=:average, weights=nothing, F1=1, F2=1, F3=1, M1=1, M2=1, M3=1)

    n = length(y)
    @assert intervals >= 1 "at least one interval per window is required"
    @assert n >= intervals + 1 "at least $(intervals + 1) nodes are required for windows of $intervals intervals"
    @assert length(T) == n - intervals "T must have the generalized displacements of the downstream node of each of the $(n - intervals) windows"
    @assert all(diff(y) .> 0) "y must be increasing"
    @assert all(size(Tk) == (6, 6) for Tk in T) "the generalized displacements of each window must be a 6 x 6 matrix"
    method in (:average, :weighted_average) || throw(ArgumentError("unknown method $method (use :average or :weighted_average)"))
    if method == :average
        weights === nothing || throw(ArgumentError("weights are only used with method = :weighted_average"))
        weights = ones(intervals)
    else
        weights = weights === nothing ? collect(1.0:intervals) : Float64.(weights)
        length(weights) == intervals || throw(ArgumentError("weights must have $intervals values (one for each position of a segment in a window), got $(length(weights))"))
        all(weights .>= 0) && sum(weights) > 0 || throw(ArgumentError("weights must be non-negative, and not all zero"))
    end

    # Generalized forces matrix
    F = diagm([F1; F2; F3; M1; M2; M3])

    # Stiffness matrix of each window
    Sw = [stiffness_matrix_from_D(Matrix{Float64}(T[w]) / F, y[w+intervals] - y[w]) for w in 1:n-intervals]

    # Weighted mean over the windows that contain each segment k, in which it is at position k - w + 1
    return map(1:n-1) do k
        ws = max(1, k - intervals + 1):min(k, n - intervals)
        sum(weights[k-w+1] * Sw[w] for w in ws) / sum(weights[k-w+1] for w in ws)
    end
end

"""
    least_squares_stiffness_matrices(; y, T, breaks=y, y_load=y[end], Cmin=nothing, λ=0, weights=nothing, F1=1, F2=1, F3=1, M1=1, M2=1, M3=1, maxiter=100_000, rtol=1e-10, details=false)

Computes the sectional stiffness matrices of the segments between the `breaks` of a beam clamped at `y[1]`, from the generalized displacements `T` of the nodes `y` under the six load cases applied at `y_load`, by a least squares fit of the compliance matrix `C = [R Z; Z' T]` of each segment, subject to `C ⪰ Cmin` (positive semidefinite), instead of the exact inversion of `stiffness_matrix` (`discretized_stiffness_matrices`) that can give stiffness matrices that are not positive definite (e.g. negative shear stiffnesses), because the transverse shear compliance of a slender beam is a tiny difference of large bending terms that is dominated by the deviations of the model from a prismatic beam

The beam is statically determinate, so by the unit load method the generalized displacements of the node at `y_j` are linear in the compliances, `T_j = ∫ P(y_j - x)' C(x) P(y_load - x) dx F` (from `y[1]` to `min(y_j, y_load)`), with `P(s) = [I 0; s ã1 I]` and the compliance constant in each segment. The fit minimizes the sum of the squared weighted residuals of all the entries of all nodes plus `λ` times the squared differences of the (scaled) compliances of neighbouring segments, which is a convex problem (unique minimum of the residual), solved by the accelerated projected gradient method (FISTA with adaptive restart), with the projection onto `C ⪰ Cmin` of each segment by its eigenvalues. With a single node at `y_load` and `breaks = [y[1], y_load]`, and no active constraint, it gives the same as `stiffness_matrix`. Returns a vector with the (6 x 6) symmetric stiffness matrix of each segment, from the root

Without `Cmin` the constraint is usually active along the transverse shear directions, and the compliance is singular (rigid shear, infinite stiffness): the data do not determine them. Set `Cmin` with the compliance of the largest admissible stiffness, e.g. `Diagonal([0, 1/GA, 1/GA, 0, 0, 0])` with `G` and `A` of the main structural member (shear stiffness at most GA), which bounds the shear stiffnesses with little or no change in the residual

# Keyword arguments
- `y` = coordinates of the nodes along the beam axis (x1-direction), increasing from the clamped node `y[1]`
- `T` = array with the (6 x 6) generalized displacements of each node, rows [u; θ] and columns [F1, F2, F3, M1, M2, M3]; the clamped node has zeros and is not used
- `breaks=y` = coordinates of the ends of the segments, increasing, from `y[1]` to `y[end]`: one stiffness matrix for each of the `length(breaks) - 1` segments. With `breaks = y` there is one for each segment between nodes, and with `breaks = [y[1], y[end]]` a single one for the whole beam, fitted to all the nodes
- `y_load=y[end]` = coordinate of the point where the loads are applied
- `Cmin=nothing` = lower bound of the compliance of every segment (6 x 6 symmetric positive semidefinite matrix, or vector with its diagonal), see above. Default: zeros
- `λ=0` = weight of the smoothness of the compliances along the beam: the squared differences between the compliances of neighbouring segments, scaled as the unknowns (see `weights`), are added to the residual with this weight
- `weights=nothing` = (6 x 6) weights of the entries of the generalized displacements of every node, row i (displacement) and column l (load case). Default `1/sqrt(|D_ii D_ll|)`, with `D` the displacements per unit load of the last node up to `y_load`, so that the residuals are relative to the direct terms. The unknowns are scaled likewise, `C = Λ C̃ Λ`, with `Λ` the square root of the diagonal of the compliance of `stiffness_matrix` at that node (at least that of `Cmin`)
- `F1=1`, `F2=1`, `F3=1`, `M1=1`, `M2=1`, `M3=1` = loads in each load case
- `maxiter=100_000` = maximum number of iterations
- `rtol=1e-10` = convergence tolerance of the relative change of the scaled compliances in an iteration
- `details=false` = return also the details of the fit: a named tuple `(S, C, T_fit, residual, iterations)` with the stiffness and compliance matrices of the segments, the generalized displacements of the nodes given by the fit (same form as `T`), the root mean square of the weighted residuals and the number of iterations
"""
function least_squares_stiffness_matrices(; y, T, breaks=y, y_load=y[end], Cmin=nothing, λ=0, weights=nothing, F1=1, F2=1, F3=1, M1=1, M2=1, M3=1, maxiter::Integer=200_000, rtol::Real=1e-10, details::Bool=false)

    n = length(y)
    nseg = length(breaks) - 1
    @assert n >= 2 "at least two nodes are required"
    @assert length(T) == n "T must have the generalized displacements of each node"
    @assert all(diff(y) .> 0) "y must be increasing"
    @assert all(size(Tk) == (6, 6) for Tk in T) "the generalized displacements of each node must be a 6 x 6 matrix"
    @assert nseg >= 1 && all(diff(breaks) .> 0) "breaks must be increasing, with at least two values"
    @assert breaks[1] ≈ y[1] && breaks[end] ≈ y[end] "breaks must go from y[1] to y[end]"
    @assert y[1] < y_load "the loads must be applied away from the clamped node"
    @assert λ >= 0 "λ must be non-negative"

    # Generalized displacements per unit load
    F = diagm([F1; F2; F3; M1; M2; M3])
    D = [Matrix{Float64}(Tk) / F for Tk in T]

    # Node used for the scaling: the last one up to the load point
    jref = findlast(<=(y_load), y)
    jref > 1 || throw(ArgumentError("no node between the clamped node and the load point"))
    Dref = D[jref]
    any(iszero, diag(Dref)) && throw(ArgumentError("the direct generalized displacements of node $jref must be non-zero"))
    W = weights === nothing ? 
        1 ./ sqrt.(abs.(diag(Dref) * diag(Dref)')) : 
        Matrix{Float64}(weights)
    size(W) == (6, 6) || throw(ArgumentError("weights must be a 6 x 6 matrix"))

    # Lower bound of the compliance
    Cmin = Cmin === nothing ? 
           zeros(6, 6) : 
           Cmin isa AbstractVector ? 
           Matrix(Diagonal(Float64.(Cmin))) : 
           Matrix{Float64}(Cmin)
    size(Cmin) == (6, 6) && issymmetric(Cmin) && eigmin(Symmetric(Cmin)) >= -eps() * max(1, opnorm(Cmin)) || throw(ArgumentError("Cmin must be a symmetric positive semidefinite 6 x 6 matrix"))

    # Scaling of the unknowns, C = Λ C̃ Λ
    C0 = diag(inv(stiffness_matrix_from_D(Dref, y[jref] - y[1])))
    Λd = sqrt.(max.(abs.(C0), diag(Cmin)))
    Λd[.!isfinite.(Λd) .| (Λd .== 0)] .= 1
    Λ = Diagonal(Λd)
    C̃min = Matrix(Symmetric(inv(Λ) * Cmin * inv(Λ)))

    # Symmetric basis of the 21 unknowns of each segment
    # TODO: the gradient step is taken in the (Euclidean) metric of these unknowns, where each off-diagonal entry counts once, whereas the projection below (by eigenvalues) is a projection in the Frobenius norm, where it counts twice. The two metrics only agree if the basis is orthonormal in the Frobenius norm, i.e. with the off-diagonal entries of the basis matrices equal to 1/sqrt(2) (and `tomat` and `tovec` scaled accordingly). With active constraints (e.g. the shear stiffness bound Cmin) the solution is then not exactly the constrained minimizer. Tested on the Technion Pazy S0 (nodes and global, nominal and refined meshes): the residual changes only from the 9th digit, the stiffness matrices of the :nodes segments by up to 0.2-2%, and the modal frequencies of S0BeamRodModal.jl not at all, so it was not adopted
    IJ = [(i, j) for j in 1:6 for i in 1:j]
    basis = map(IJ) do (i, j)
        B = zeros(6, 6); B[i, j] = B[j, i] = 1; Λ * B * Λ
    end
    tomat(c) = (M = zeros(6, 6); for (k, (i, j)) in enumerate(IJ); M[i, j] = M[j, i] = c[k]; end; M)
    tovec(M) = [M[i, j] for (i, j) in IJ]

    # Generalized forces transmitted by the section at distance s from the point of the load (or of the unit load)
    a1tilde = [0 0 0; 0 0 -1; 0 1 0]
    P(s) = [I(3) zeros(3, 3); s*a1tilde I(3)]

    # Linear operator of the weighted generalized displacements of the nodes (all but the clamped one) on the scaled unknowns: the integrand is quadratic in x, so the 2-point Gauss rule is exact
    nodes = 2:n
    A = zeros(36 * length(nodes), 21 * nseg)
    b = zeros(36 * length(nodes))
    for (r, j) in enumerate(nodes)
        rows = 36(r-1)+1:36r
        b[rows] = vec(W .* D[j])
        for s in 1:nseg
            x0, x1 = breaks[s], min(breaks[s+1], y[j], y_load)
            x1 > x0 || continue
            xg = (x0 + x1) / 2 .+ (x1 - x0) / 2 .* (-1, 1) ./ sqrt(3)
            for k in 1:21
                K = (x1 - x0) / 2 * sum(P(y[j] - x)' * basis[k] * P(y_load - x) for x in xg)
                A[rows, 21(s-1)+k] = vec(W .* K)
            end
        end
    end

    # Smoothness: differences between the unknowns of neighbouring segments
    Δ = zeros(21 * (nseg - 1), 21 * nseg)
    for s in 1:nseg-1, k in 1:21
        Δ[21(s-1)+k, 21(s-1)+k] = -1
        Δ[21(s-1)+k, 21s+k] = 1
    end

    # Projection of the unknowns of each segment onto C̃ ⪰ C̃min
    function project(v)
        v = copy(v)
        for s in 1:nseg
            idx = 21(s-1)+1:21s
            ev = eigen(Symmetric(tomat(v[idx]) - C̃min))
            v[idx] = tovec(ev.vectors * Diagonal(max.(ev.values, 0)) * ev.vectors' + C̃min)
        end
        return v
    end

    # FISTA with adaptive (gradient) restart, from the projection of the unconstrained solution
    H = A' * A + λ * Δ' * Δ
    g = A' * b
    Lip = opnorm(H)
    c = project(pinv(H) * g)
    z, t = copy(c), 1.0
    iterations = maxiter
    for it in 1:maxiter
        cnew = project(z - (H * z - g) / Lip)
        if dot(z - cnew, cnew - c) > 0
            z, t = copy(cnew), 1.0
        else
            tnew = (1 + sqrt(1 + 4t^2)) / 2
            z = cnew + (t - 1) / tnew * (cnew - c)
            t = tnew
        end
        converged = norm(cnew - c) <= rtol * max(norm(cnew), 1)
        c = cnew
        if converged
            iterations = it
            break
        end
    end
    iterations == maxiter && @warn "least_squares_stiffness_matrices did not converge in $maxiter iterations"

    # Compliance and stiffness matrices of the segments
    C = [Matrix(Symmetric(Λ * tomat(c[21(s-1)+1:21s]) * Λ)) for s in 1:nseg]
    for s in 1:nseg
        e = eigvals(Symmetric(tomat(c[21(s-1)+1:21s])))
        e[1] <= 1e-10 * e[end] && @warn "the compliance matrix of segment $s is singular (smallest scaled eigenvalue $(e[1])): its stiffness matrix is not finite. Set Cmin to bound it"
    end
    S = [Matrix(Symmetric(inv(Symmetric(Ck)))) for Ck in C]

    details || return S

    Afit = A * c
    T_fit = [zeros(6, 6) for _ in 1:n]
    for (r, j) in enumerate(nodes)
        T_fit[j] = reshape(Afit[36(r-1)+1:36r], 6, 6) ./ W * F
    end
    residual = sqrt(sum(abs2, Afit - b) / length(b))
    return (S=S, C=C, T_fit=T_fit, residual=residual, iterations=iterations)
end

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


"""
    stiffness_matrix_normalized_derivatives(; kwargs...)

Computes the normalized Jacobian matrix of the sectional stiffness matrix with respect to the generalized displacements resulting from the loads

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
function stiffness_matrix_normalized_derivatives(; F1=1,F2=1,F3=1,M1=1,M2=1,M3=1,L,u_F1,u_F2,u_F3,u_M1,u_M2,u_M3,θ_F1,θ_F2,θ_F3,θ_M1,θ_M2,θ_M3)

    # Array of generalized displacements
    genDisp = vcat(u_F1,θ_F1,u_F2,θ_F2,u_F3,θ_F3,u_M1,θ_M1,u_M2,θ_M2,u_M3,θ_M3)

    # Jacobian matrix
    J = ForwardDiff.jacobian(x -> stiffness_matrix_wrapper(x,F1,F2,F3,M1,M2,M3,L), genDisp)

    # Stiffness
    S_vec = stiffness_matrix_wrapper(genDisp, F1, F2, F3, M1, M2, M3, L)

    # Normalization
    J_norm = zeros(eltype(J), size(J))
    for i in axes(J, 1)
        for j in axes(J, 2)
            if S_vec[i] != 0.0
                J_norm[i, j] = J[i, j] * genDisp[j] / S_vec[i]
            else
                J_norm[i, j] = 0.0
            end
        end
    end

    return J_norm, J
end

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
    D11 = L*R + L^2/2*Z*a1tilde - L^2/2*a1tilde*Z' - L^3/3*a1tilde*T*a1tilde
    D12 = L*Z - L^2/2*a1tilde*T
    D21 = L*Z' + L^2/2*T*a1tilde
    D22 = L*T

    D = [D11 D12; D21 D22]

    # Generalized displacements
    uθ = D * F

    return uθ
end