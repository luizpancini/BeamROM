# ------------------------------------------------------------------------------------------------------------
# Pazy wing (Technion model)
# ------------------------------------------------------------------------------------------------------------

# Properties of the CQUAD4 elements of the ribs of the Pazy wing
const PAZY_RIB_PIDS = (4, 10, 12, 200004, 200010)

"""
    rib_stations(gfem; rib_pids=PAZY_RIB_PIDS, tol=1e-5) -> Vector{Float64}

Spanwise (y) positions of the ribs of the Pazy wing, sorted from the root: the planes (constant y) of the shell elements with properties `rib_pids`. Elements not contained in a plane of constant y (e.g. the root transition elements) are not ribs. Planes closer than `tol` are the same rib.
"""
function rib_stations(gfem::GFEM; rib_pids=PAZY_RIB_PIDS, tol::Real=1e-5)
    ys = Float64[]
    for (_, pid, gids) in gfem.shells
        pid in rib_pids || continue
        y = [gfem.grids[g][2] for g in gids]
        maximum(y) - minimum(y) < tol && push!(ys, sum(y) / length(y))
    end
    isempty(ys) && throw(ArgumentError("no planar rib elements with properties $(collect(rib_pids)) in the model"))
    sort!(ys)
    ribs = [ys[1]]
    for y in ys[2:end]
        y - ribs[end] > tol && push!(ribs, y)
    end
    return ribs
end

"""
    spar_reference_axis(gfem, nodes; reference_axis, sweep, spar_pid=100003, rib_pids=PAZY_RIB_PIDS, skip_root=2, skip_tip=3) -> Vector{Vector{Float64}}

Beam nodes along the reference axis of the spar of the Pazy wing: the spar is modeled by the CQUAD4 elements of property `spar_pid` (the 60 × 2.25 mm aluminum plate, on the plane z = 0), so the axis lies on its mid-surface (z at the mid-thickness plane). The wing is along y, with the chord along x. The axis goes from the root to the tip of the spar. `nodes` is either the number of equally spaced nodes, a vector of fractions of the span (as in `straight_axis`), `:ribs` or `:ribs_halved`. With `:ribs` the nodes are the root of the spar, the ribs (see `rib_stations`) except the first `skip_root` ones, close to the root and lumped in the root node, and the last `skip_tip` ones, close to the tip and lumped in the tip node, and the tip of the spar. `:ribs_halved` adds a node at the middle of each segment between consecutive `:ribs` nodes, except the last one (the short segment to the tip of the spar, with the 4 ribs close to the tip), so it has twice as many segments minus one; each of these nodes is moved to the nearest row of spar grids that spans the chord of the spar (the lower one in a tie; rows with one or a few grids, as near the root of the S10 model, are skipped), because the stiffness workflows need the section of the spar at the nodes. The chordwise (x) position of the axis is the fraction `reference_axis` (required, between 0 and 1) of the chord of the wing, measured from the minimum x of the shell grids, with the chord being the x extent of the shell grids (e.g. 0.25 is the quarter-chord if x grows toward the trailing edge). `sweep` is the sweep angle [deg] of a swept wing whose sections (planes of constant y) are shifted along x by `y tan(sweep)` (e.g. the S10 and S20 models, with the spar sheared along the axis; 0 for the straight wing, S0; required): the axis is then inclined by `sweep` from the y axis, and the chordwise position (at `reference_axis` of the chord) is taken in the unswept coordinates `x - y tan(sweep)`, i.e. at a constant fraction of the chord along the span.
"""
function spar_reference_axis(gfem::GFEM, nodes::Union{Integer,AbstractVector{<:Real},Symbol}; reference_axis::Real, spar_pid::Integer=100003, rib_pids=PAZY_RIB_PIDS, skip_root::Integer=2, skip_tip::Integer=3, sweep::Real)
    0 <= reference_axis <= 1 || throw(ArgumentError("reference_axis must be a fraction of the chord between 0 and 1, got $reference_axis"))
    gids = unique(vcat([g for (_, pid, g) in gfem.shells if pid == spar_pid]...))
    isempty(gids) && throw(ArgumentError("no shell elements with property $spar_pid in the model"))
    P = reduce(hcat, gfem.grids[g] for g in gids)
    # Chordwise coordinate of the sections without the sweep
    tanΛ = tand(sweep)
    unswept(p) = p[1] - tanΛ * p[2]
    xs = [unswept(gfem.grids[g]) for (_, _, gids) in gfem.shells for g in gids]
    xc = minimum(xs) + reference_axis * (maximum(xs) - minimum(xs))
    zc = (minimum(P[3, :]) + maximum(P[3, :])) / 2
    y0, y1 = minimum(P[2, :]), maximum(P[2, :])
    if nodes isa Symbol
        halved = nodes == :ribs_halved
        nodes in (:ribs, :ribs_halved) || throw(ArgumentError("unknown node distribution $nodes (use :ribs or :ribs_halved)"))
        ribs = rib_stations(gfem; rib_pids)
        interior = ribs[skip_root+1:end-skip_tip]
        all(y -> y0 < y < y1, interior) || throw(ArgumentError("the selected ribs are not inside the spar"))
        nodes = ([y0; interior; y1] .- y0) / (y1 - y0)
        if halved
            # Each midpoint goes to the nearest row of spar grids (the lower one in a tie), so that the node is at a section of the spar. Only the rows that span the chord of the spar are used: the irregular mesh near the root has rows with one or a few grids (e.g. in the S10 model), at which the rigid body motion of the section is not determined
            rows = sort(P[2, :])
            rows = rows[[true; diff(rows) .> 1e-5]]
            ux = [unswept(P[:, k]) for k in axes(P, 2)]
            sparChord = maximum(ux) - minimum(ux)
            spans(y) = (u = ux[abs.(P[2, :] .- y) .< 1e-5]; maximum(u) - minimum(u) >= 0.9 * sparChord)
            rows = filter(spans, rows)
            mids = (nodes[1:end-2] + nodes[2:end-1]) / 2
            snapped = map(mids) do m
                ym = y0 + m * (y1 - y0)
                (rows[argmin(round.(abs.(rows .- ym), digits=9))] - y0) / (y1 - y0)
            end
            all(k -> nodes[k] < snapped[k] < nodes[k+1], eachindex(snapped)) || throw(ArgumentError("a halved segment has no row of spar grids between its nodes"))
            nodes = sort!([nodes; snapped])
        end
    end
    return straight_axis([xc + tanΛ * y0, y0, zc], [xc + tanΛ * y1, y1, zc], nodes)
end

"""
    pazy_beam_mass_model(bdf; reference_axis, sweep, nodes=16, spar_pid=100003, rib_pids=PAZY_RIB_PIDS, skip_root=2, skip_tip=3, kwargs...)

Equivalent beam mass model of the Pazy wing Nastran model `bdf` (e.g. `S0.bdf`, with shell, CBEAM and CONM2 elements), condensed to the reference axis of the aluminum spar (see `spar_reference_axis` for the node distribution `nodes`, which can follow the ribs with `nodes=:ribs` or `:ribs_halved`). The wing is along y and the chord along x, so the grids are assigned to the nearest beam node along y, and the offsets and inertia tensors are reported in the Nastran basic frame (x chordwise, y spanwise, z normal): Jyy is the torsional inertia, Jxx and Jzz the inertias about the chord and the normal axes. The root grids belong to the first node and the grids beyond the spar tip (tip rods and tip mass) to the last one. `reference_axis` (required, between 0 and 1) places the reference axis at that fraction of the chord. `sweep` is the sweep angle [deg] of a swept wing with the sections in planes of constant y (e.g. the S10 and S20 models, 0 for the straight wing, S0; required, see `spar_reference_axis`): the reference axis is inclined by that angle from the y axis, but the grids are still assigned along y and the offsets and tensors are still in the Nastran basic frame, so Jyy is the inertia about the y axis, not about the swept axis. The keywords `kwargs` are those of `beam_mass_model`, except `span_direction`, which is the y axis.
"""
function pazy_beam_mass_model(bdf::AbstractString; reference_axis::Real, nodes::Union{Integer,AbstractVector{<:Real},Symbol}=16, spar_pid::Integer=100003, rib_pids=PAZY_RIB_PIDS, skip_root::Integer=2, skip_tip::Integer=3, sweep::Real, kwargs...)
    xref = spar_reference_axis(read_gfem(bdf), nodes; spar_pid, rib_pids, skip_root, skip_tip, reference_axis, sweep)
    return beam_mass_model(bdf, xref; kwargs..., span_direction=[0.0, 1.0, 0.0])
end
