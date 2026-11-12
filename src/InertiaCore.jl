# Generation of the equivalent beam mass model following the methodology of the "Equivalent Beam Mass Model":
#
#   The GFEM is represented as a collection of lumped masses m_i located at x_i + δ_i (x_i = GFEM grid position, δ_i = mass offset) with lumped inertia tensors J_i about x_i + δ_i, all in the Nastran frame. The equivalent mass, offset and inertia tensor at each beam reference-axis node are obtained by  summing the contributions of the GFEM grids in its nearest neighborhood:
#
#   m_ref = Σ m_i,   
#   δ_ref = (1/m_ref) Σ m_i (x_i + δ_i − x_ref),   
#   J_ref = Σ (J_i + J_i,tr)
#
# 1. Lumped GFEM (equivalent to the diagonal blocks of Nastran's lumped mass matrix, PARAM,COUPMASS,-1):
#    - CHEXA (8 nodes, PSOLID/MAT1): mass ρ V distributed to the grids by the row sum of the consistent mass matrix, ρ ∫ Nᵢ dV (2×2×2 Gauss quadrature of the trilinear isoparametric map; V/8 per grid for parallelepipeds), no rotational inertia;
#    - CQUAD4/CTRIA3: mass (ρ t + NSM) A distributed to the grids by the row sum of the consistent mass matrix, ∫ Nᵢ dA (equal split for parallelograms/triangles), no rotational inertia;
#    - CBEAM/CBAR: half of ∫(ρ A + NSM) dx at each end grid, placed at the PBEAM neutral axis offset (N1, N2) (PBAR: on the element axis); the torsional rotary inertia ρ (I1 + I2) L / 2 about the element axis is also lumped at each end;
#    - CONM2: mass, offset and inertia tensor about the center of mass taken directly from the card (Nastran's I21, I31, I32 are products of inertia, so the tensor off-diagonals are their negatives).
# 2. Reference axis: specified by the user as the list of beam node coordinates (e.g. `straight_axis`). For simple prismatic models it is the elastic axis, which coincides with the shear center.
# 3. Nearest neighborhood: each GFEM grid is assigned to the beam node with the closest coordinate along the span direction (by default the direction of the reference axis, but it can be set, e.g. to the unswept span of a swept wing), i.e. the structure is sliced into bins centered at the beam nodes (the Euclidean distance to a swept axis gives a different and wrong split between the last nodes). The mass of the grids that lie exactly midway between beam nodes (up to round-off) is split equally between them.
# 4. The inertia tensors are reported about the center of mass of each beam node (x_ref + δ_ref), in the NASTRAN frame and with the tensor sign convention (J_xy = −∫ x y dm).

using LinearAlgebra
using Printf

# ------------------------------------------------------------------------------------------------------------
# Nastran bulk data parsing
# ------------------------------------------------------------------------------------------------------------

"""
    nastran_real(s) -> Float64

Parses a Nastran real field, including the shorthand exponent notation without `E` (e.g. `4.3444001-5`, `1.-8`, `1.4+9`).
"""
function nastran_real(s::AbstractString)
    t = replace(strip(s), 'D' => 'E', 'd' => 'e')
    isempty(t) && throw(ArgumentError("blank real field"))
    if !occursin(r"[eE]", t)
        m = match(r"^([+-]?[0-9]*\.?[0-9]*)([+-][0-9]+)$", t)
        if m !== nothing
            t = m.captures[1] * "e" * m.captures[2]
        end
    end
    # Julia does not accept a trailing dot before the exponent (e.g. "6.e-5")
    t = replace(t, r"\.([eE])" => s".0\1")
    return parse(Float64, t)
end

real_or(fields, i, default) = (i <= length(fields) && !isempty(strip(fields[i]))) ? nastran_real(fields[i]) : default
int_or(fields, i, default) = (i <= length(fields) && !isempty(strip(fields[i]))) ? parse(Int, strip(fields[i])) : default
str_or(fields, i, default) = (i <= length(fields) && !isempty(strip(fields[i]))) ? String(strip(fields[i])) : default

substr(line, a, b) = a > length(line) ? "" : line[a:min(b, length(line))]

"""
    split_bulk_line(line) -> (name, data, iscontinuation)

Splits one physical bulk-data line into its name and data fields, handling small-field (8 columns), large-field (16 columns, `*` in the first field) and free-field (comma-separated) formats.
"""
function split_bulk_line(line::AbstractString)
    if occursin(',', line)
        parts = strip.(split(line, ','))
        name = parts[1]
        data = String.(parts[2:end])
    else
        name = strip(substr(line, 1, 8))
        if occursin('*', name)
            data = [String(strip(substr(line, 9 + 16(k - 1), 8 + 16k))) for k in 1:4]
        else
            data = [String(strip(substr(line, 9 + 8(k - 1), 8 + 8k))) for k in 1:8]
        end
    end
    iscont = isempty(name) || startswith(name, '+') || startswith(name, '*')
    return String(replace(name, "*" => "")), data, iscont
end

"""
    read_bulk_lines(path) -> Vector{String}

Reads the bulk data lines of a Nastran input file, resolving `INCLUDE` statements (relative to the file folder),
skipping comments and, if present, everything outside the `BEGIN BULK`/`ENDDATA` section.
"""
function read_bulk_lines(path::AbstractString)
    raw = readlines(path)
    ibegin = findfirst(l -> occursin(r"^\s*BEGIN\s+BULK"i, l), raw)
    iend = findfirst(l -> occursin(r"^\s*ENDDATA"i, l), raw)
    first = ibegin === nothing ? 1 : ibegin + 1
    last = iend === nothing ? length(raw) : iend - 1
    lines = String[]
    for l in raw[first:last]
        s = rstrip(l)
        (isempty(strip(s)) || startswith(lstrip(s), '$')) && continue
        m = match(r"^\s*INCLUDE\s+'?([^']+?)'?\s*$"i, s)
        if m !== nothing
            append!(lines, read_bulk_lines(joinpath(dirname(path), strip(m.captures[1]))))
        else
            push!(lines, s)
        end
    end
    return lines
end

"""
    read_cards(path) -> Vector{Pair{String,Vector{String}}}

Assembles the logical cards (name => data fields) of a Nastran input file, concatenating continuation lines.
"""
function read_cards(path::AbstractString)
    cards = Pair{String,Vector{String}}[]
    for line in read_bulk_lines(path)
        name, data, iscont = split_bulk_line(line)
        if iscont
            isempty(cards) && error("continuation line without parent card: $line")
            append!(cards[end].second, data)
        else
            push!(cards, uppercase(name) => data)
        end
    end
    return cards
end

# ------------------------------------------------------------------------------------------------------------
# GFEM data structures
# ------------------------------------------------------------------------------------------------------------

struct PBeam
    mid::Int
    xstations::Vector{Float64}   # x/xb of the property stations
    A::Vector{Float64}           # area at the stations
    Ip::Vector{Float64}          # polar area moment of inertia I1 + I2 at the stations
    nsm::Vector{Float64}         # nonstructural mass per length at the stations
    m1::NTuple{2,Float64}        # NSM center of gravity offset (y,z) at end A
    m2::NTuple{2,Float64}        # ... at end B
    n1::NTuple{2,Float64}        # neutral axis offset (y,z) at end A
    n2::NTuple{2,Float64}        # ... at end B
end

struct GFEM
    grids::Dict{Int,Vector{Float64}}
    shells::Vector{Tuple{Int,Int,Vector{Int}}} # (eid, pid, grids)
    solids::Vector{Tuple{Int,Int,Vector{Int}}} # (eid, pid, grids) of the CHEXA elements
    beams::Vector{Tuple{Int,Int,Int,Int,Union{Int,Vector{Float64}}}} # (eid, pid, ga, gb, orientation) of the CBEAM elements
    bars::Vector{Tuple{Int,Int,Int,Int,Union{Int,Vector{Float64}}}} # same as beams, for the CBAR elements
    conm2::Vector{Tuple{Int,Int,Int,Float64,Vector{Float64},Matrix{Float64}}}  # (eid, grid, cid, m, offset, J)
    pshell::Dict{Int,NamedTuple{(:mid1, :mid2, :t, :nsm),Tuple{Int,Int,Float64,Float64}}}
    psolid::Dict{Int,Int}        # pid => mid
    pbeam::Dict{Int,PBeam}
    pbar::Dict{Int,PBeam}        # PBAR cards, represented as constant-section PBeam properties without offsets
    rho::Dict{Int,Float64}
end

# Mass-carrying cards that are not handled (their mass would be silently ignored)
const UNSUPPORTED_MASS_CARDS = ("CTETRA", "CPENTA", "CPYRAM", "CROD", "CONROD", "CTUBE", "CMASS1", "CMASS2", "CMASS3", "CMASS4", "CONM1", "CQUAD8", "CTRIA6", "CQUADR", "CTRIAR", "CBEAM3", "CBARAO", "PBARL", "PBEAML")

"""
    parse_pbeam(f) -> PBeam

Parses the fields of a PBEAM card, with an arbitrary number of intermediate stations (SO = YES, YESA or NO).
"""
function parse_pbeam(f::Vector{String})
    mid = int_or(f, 2, 0)
    xs = [0.0]
    A = [real_or(f, 3, 0.0)]
    I1 = [real_or(f, 4, 0.0)]
    I2 = [real_or(f, 5, 0.0)]
    nsm = [real_or(f, 8, 0.0)]
    i = 17
    # Blank section properties at intermediate stations/end B default to the end A values
    while i <= length(f) && uppercase(strip(f[i])) in ("YES", "YESA", "NO")
        so = uppercase(strip(f[i]))
        push!(xs, real_or(f, i + 1, 1.0))
        push!(A, real_or(f, i + 2, A[1]))
        push!(I1, real_or(f, i + 3, I1[1]))
        push!(I2, real_or(f, i + 4, I2[1]))
        push!(nsm, real_or(f, i + 7, nsm[1]))
        i += so == "YES" ? 16 : 8
    end
    # Trailing lines: K1 K2 S1 S2 NSI(A) NSI(B) CW(A) CW(B) / M1(A) M2(A) M1(B) M2(B) N1(A) N2(A) N1(B) N2(B)
    m1 = (real_or(f, i + 8, 0.0), real_or(f, i + 9, 0.0))
    m2 = (real_or(f, i + 10, m1[1]), real_or(f, i + 11, m1[2]))
    n1 = (real_or(f, i + 12, 0.0), real_or(f, i + 13, 0.0))
    n2 = (real_or(f, i + 14, n1[1]), real_or(f, i + 15, n1[2]))
    # If the beam has a single station, Nastran uses the end A values at end B
    if length(xs) == 1
        push!(xs, 1.0); push!(A, A[1]); push!(I1, I1[1]); push!(I2, I2[1]); push!(nsm, nsm[1])
    end
    return PBeam(mid, xs, A, I1 .+ I2, nsm, m1, m2, n1, n2)
end

"""
    parse_pbar(f) -> PBeam

Parses the fields of a PBAR card (constant section, no neutral axis offset). The polar inertia ρ(I1 + I2) is used for the torsional rotary inertia.
"""
function parse_pbar(f::Vector{String})
    A, I1, I2, nsm = real_or(f, 3, 0.0), real_or(f, 4, 0.0), real_or(f, 5, 0.0), real_or(f, 7, 0.0)
    return PBeam(int_or(f, 2, 0), [0.0, 1.0], [A, A], [I1 + I2, I1 + I2], [nsm, nsm], (0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0))
end

# Orientation of a CBEAM/CBAR: vector (decimal point) or G0 grid, and a check for the unsupported end offsets
function line_element(f::Vector{String}, kind::AbstractString)
    orient = occursin('.', str_or(f, 5, "0.")) ? [real_or(f, k, 0.0) for k in 5:7] : int_or(f, 5, 0)
    any(k -> real_or(f, k, 0.0) != 0.0, 11:16) && error("$kind $(f[1]) has end offsets (not supported)")
    return (int_or(f, 1, 0), int_or(f, 2, 0), int_or(f, 3, 0), int_or(f, 4, 0), orient)
end

"""
    read_gfem(path) -> GFEM

Reads the GRID, CHEXA, CQUAD4, CTRIA3, CBEAM, CBAR, CONM2, PSOLID, PSHELL, PBEAM, PBAR and MAT1 cards of a Nastran input file. A warning is issued for mass-carrying cards that are not supported.
"""
function read_gfem(path::AbstractString)
    gfem = GFEM(Dict(), [], [], [], [], [], Dict(), Dict(), Dict(), Dict(), Dict())
    ignored = Dict{String,Int}()
    for (name, f) in read_cards(path)
        if name == "GRID"
            int_or(f, 2, 0) == 0 || error("GRID $(f[1]) uses a non-basic coordinate system (not supported)")
            gfem.grids[int_or(f, 1, 0)] = [real_or(f, 3, 0.0), real_or(f, 4, 0.0), real_or(f, 5, 0.0)]
        elseif name == "CQUAD4"
            push!(gfem.shells, (int_or(f, 1, 0), int_or(f, 2, 0), [int_or(f, k, 0) for k in 3:6]))
        elseif name == "CTRIA3"
            push!(gfem.shells, (int_or(f, 1, 0), int_or(f, 2, 0), [int_or(f, k, 0) for k in 3:5]))
        elseif name == "CHEXA"
            # Only the 8 corner grids are supported (a 20-node CHEXA has midside grids in fields 11-20)
            any(k -> int_or(f, k, 0) != 0, 11:length(f)) && error("CHEXA $(f[1]) has midside grids (not supported)")
            push!(gfem.solids, (int_or(f, 1, 0), int_or(f, 2, 0), [int_or(f, k, 0) for k in 3:10]))
        elseif name == "CBEAM"
            push!(gfem.beams, line_element(f, "CBEAM"))
        elseif name == "CBAR"
            push!(gfem.bars, line_element(f, "CBAR"))
        elseif name == "CONM2"
            cid = int_or(f, 3, 0)
            cid in (0, -1) || error("CONM2 $(f[1]) uses a non-basic coordinate system (not supported)")
            m = real_or(f, 4, 0.0)
            offset = [real_or(f, k, 0.0) for k in 5:7]
            I11, I21, I22, I31, I32, I33 = (real_or(f, k, 0.0) for k in 9:14)
            # Nastran's I21, I31, I32 are mass products of inertia: the tensor off-diagonals carry a minus sign
            J = [I11 -I21 -I31; -I21 I22 -I32; -I31 -I32 I33]
            push!(gfem.conm2, (int_or(f, 1, 0), int_or(f, 2, 0), cid, m, offset, J))
        elseif name == "PSHELL"
            gfem.pshell[int_or(f, 1, 0)] = (mid1=int_or(f, 2, 0), mid2=int_or(f, 4, 0), t=real_or(f, 3, 0.0), nsm=real_or(f, 8, 0.0))
        elseif name == "PSOLID"
            gfem.psolid[int_or(f, 1, 0)] = int_or(f, 2, 0)
        elseif name == "PBEAM"
            gfem.pbeam[int_or(f, 1, 0)] = parse_pbeam(f)
        elseif name == "PBAR"
            gfem.pbar[int_or(f, 1, 0)] = parse_pbar(f)
        elseif name == "MAT1"
            gfem.rho[int_or(f, 1, 0)] = real_or(f, 5, 0.0)
        elseif name in UNSUPPORTED_MASS_CARDS
            ignored[name] = get(ignored, name, 0) + 1
        end
    end
    for (name, n) in ignored
        @warn "$n $name card(s) ignored: their mass is not included in the GFEM mass model"
    end
    return gfem
end

# ------------------------------------------------------------------------------------------------------------
# Lumped mass representation of the GFEM
# ------------------------------------------------------------------------------------------------------------

"""
    LumpedMass

Lumped mass `m` located at `p` (basic frame) with inertia tensor `J` about `p`, attached to GFEM grid `grid`.
"""
struct LumpedMass
    grid::Int
    m::Float64
    p::Vector{Float64}
    J::Matrix{Float64}
end

tilde(v) = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]

# Inertia tensor transport term of mass m from its center of mass to a point displaced by r
transport(m, r) = m * (dot(r, r) * I(3) - r * r')

triangle_area(a, b, c) = 0.5 * norm(cross(b - a, c - a))

"""
    shell_nodal_areas(coords, lumping) -> Vector{Float64}

Tributary areas of the grids of a 3- or 4-node shell element, whose sum is the element area. For triangles all lumping schemes give A/3 per grid. For quadrilaterals (bilinear isoparametric map, 2×2 Gauss quadrature):
- `:equal`  → A/4 per grid;
- `:rowsum` → ∫ Nᵢ dA (row sum of the consistent mass matrix);
- `:hrz`    → A ∫ Nᵢ² dA / Σⱼ ∫ Nⱼ² dA (diagonal scaling of the consistent mass matrix).
"""
function shell_nodal_areas(coords::Vector{Vector{Float64}}, lumping::Symbol)
    if length(coords) == 3
        return fill(triangle_area(coords...) / 3, 3)
    end
    g = 1 / sqrt(3)
    N(ξ, η) = [(1 - ξ) * (1 - η), (1 + ξ) * (1 - η), (1 + ξ) * (1 + η), (1 - ξ) * (1 + η)] / 4
    dNξ(η) = [-(1 - η), (1 - η), (1 + η), -(1 + η)] / 4
    dNη(ξ) = [-(1 - ξ), -(1 + ξ), (1 + ξ), (1 - ξ)] / 4
    intN = zeros(4)
    intN2 = zeros(4)
    for ξ in (-g, g), η in (-g, g)
        xξ = sum(dNξ(η)[i] * coords[i] for i in 1:4)
        xη = sum(dNη(ξ)[i] * coords[i] for i in 1:4)
        detJ = norm(cross(xξ, xη))
        n = N(ξ, η)
        intN .+= n * detJ
        intN2 .+= n .^ 2 * detJ
    end
    A = sum(intN)
    if lumping == :equal
        return fill(A / 4, 4)
    elseif lumping == :rowsum
        return intN
    elseif lumping == :hrz
        return A * intN2 / sum(intN2)
    else
        throw(ArgumentError("unknown shell lumping scheme $lumping"))
    end
end

"""
    solid_nodal_volumes(coords, lumping) -> Vector{Float64}

Tributary volumes of the 8 grids of a CHEXA element (trilinear isoparametric map, 2×2×2 Gauss quadrature), whose sum is the element volume. Grids in Nastran order (1-4 bottom face, 5-8 top face, counter-clockwise):
- `:equal`  → V/8 per grid;
- `:rowsum` → ∫ Nᵢ dV (row sum of the consistent mass matrix);
- `:hrz`    → V ∫ Nᵢ² dV / Σⱼ ∫ Nⱼ² dV (diagonal scaling of the consistent mass matrix).
"""
function solid_nodal_volumes(coords::Vector{Vector{Float64}}, lumping::Symbol)
    length(coords) == 8 || throw(ArgumentError("only 8-node hexahedra are supported"))
    ξi = (-1, 1, 1, -1, -1, 1, 1, -1)
    ηi = (-1, -1, 1, 1, -1, -1, 1, 1)
    ζi = (-1, -1, -1, -1, 1, 1, 1, 1)
    g = 1 / sqrt(3)
    intN = zeros(8)
    intN2 = zeros(8)
    for ξ in (-g, g), η in (-g, g), ζ in (-g, g)
        n = [(1 + ξ * ξi[i]) * (1 + η * ηi[i]) * (1 + ζ * ζi[i]) / 8 for i in 1:8]
        dξ = [ξi[i] * (1 + η * ηi[i]) * (1 + ζ * ζi[i]) / 8 for i in 1:8]
        dη = [(1 + ξ * ξi[i]) * ηi[i] * (1 + ζ * ζi[i]) / 8 for i in 1:8]
        dζ = [(1 + ξ * ξi[i]) * (1 + η * ηi[i]) * ζi[i] / 8 for i in 1:8]
        Jac = hcat(sum(dξ[i] * coords[i] for i in 1:8), sum(dη[i] * coords[i] for i in 1:8), sum(dζ[i] * coords[i] for i in 1:8))
        detJ = abs(det(Jac))
        intN .+= n * detJ
        intN2 .+= n .^ 2 * detJ
    end
    V = sum(intN)
    if lumping == :equal
        return fill(V / 8, 8)
    elseif lumping == :rowsum
        return intN
    elseif lumping == :hrz
        return V * intN2 / sum(intN2)
    else
        throw(ArgumentError("unknown solid lumping scheme $lumping"))
    end
end

"""
    beam_frame(xa, xb, orient, grids) -> (x̂, ŷ, ẑ, L)

Element coordinate system of a CBEAM/CBAR from its end coordinates and orientation vector (or G0 grid).
"""
function beam_frame(xa, xb, orient, grids)
    L = norm(xb - xa)
    x̂ = (xb - xa) / L
    v = orient isa Int ? grids[orient] - xa : orient
    ŷ = v - dot(v, x̂) * x̂
    ŷ /= norm(ŷ)
    return x̂, ŷ, cross(x̂, ŷ), L
end

"""
    lump_gfem(gfem; shell_lumping=:rowsum, solid_lumping=:rowsum, beam_axis_offsets=true, beam_torsional_inertia=true) -> Vector{LumpedMass}

Builds the lumped mass representation of the GFEM, equivalent to Nastran's lumped (PARAM,COUPMASS,-1) mass matrix:
- CHEXA elements: mass ρ V distributed to the element grids according to `solid_lumping` (see `solid_nodal_volumes`), no rotational inertia;
- shell elements: mass (ρ t + NSM) A distributed to the element grids according to `shell_lumping` (see `shell_nodal_areas`), no rotational inertia;
- beam elements (CBEAM and CBAR): mass ∫(ρ A + NSM) dx lumped at the end grids (consistent with linear interpolation along the element). If `beam_axis_offsets` is true the structural and nonstructural masses are placed at the neutral axis (N1, N2) and NSM center of gravity (M1, M2) offsets of the PBEAM (a PBAR has none). If `beam_torsional_inertia` is true, the rotary inertia about the element axis ∫ρ(I1 + I2)dx is also lumped at the end grids (no bending rotary inertia), as MSC Nastran does for the CBEAM;
- CONM2: mass, offset and inertia tensor taken directly from the card.
"""
function lump_gfem(gfem::GFEM; shell_lumping::Symbol=:rowsum, solid_lumping::Symbol=:rowsum, beam_axis_offsets::Bool=true, beam_torsional_inertia::Bool=true)
    lumped = LumpedMass[]
    Z = zeros(3, 3)
    for (eid, pid, gids) in gfem.solids
        ρ = get(gfem.rho, gfem.psolid[pid], 0.0)
        ρ == 0.0 && continue
        coords = [gfem.grids[g] for g in gids]
        volumes = solid_nodal_volumes(coords, solid_lumping)
        for (g, x, v) in zip(gids, coords, volumes)
            push!(lumped, LumpedMass(g, ρ * v, x, Z))
        end
    end
    for (eid, pid, gids) in gfem.shells
        prop = gfem.pshell[pid]
        mid = prop.mid1 != 0 ? prop.mid1 : prop.mid2
        μ = get(gfem.rho, mid, 0.0) * prop.t + prop.nsm
        μ == 0.0 && continue
        coords = [gfem.grids[g] for g in gids]
        areas = shell_nodal_areas(coords, shell_lumping)
        for (g, x, a) in zip(gids, coords, areas)
            push!(lumped, LumpedMass(g, μ * a, x, Z))
        end
    end
    for (elements, properties) in ((gfem.beams, gfem.pbeam), (gfem.bars, gfem.pbar)), (eid, pid, ga, gb, orient) in elements
        prop = properties[pid]
        ρ = get(gfem.rho, prop.mid, 0.0)
        xa, xb = gfem.grids[ga], gfem.grids[gb]
        x̂, ŷ, ẑ, L = beam_frame(xa, xb, orient, gfem.grids)
        # Structural and nonstructural mass per length, and torsional mass inertia per length, at the stations
        μs = ρ .* prop.A
        μn = prop.nsm
        μt = ρ .* prop.Ip
        # Lumping with linear shape functions: m_A = L∫μ(1-ξ)dξ, m_B = L∫μ ξ dξ (trapezoidal between stations)
        function lump_ends(μ)
            mA = mB = 0.0
            for k in 1:length(prop.xstations)-1
                ξ1, ξ2 = prop.xstations[k], prop.xstations[k+1]
                μ1, μ2 = μ[k], μ[k+1]
                Δ = ξ2 - ξ1
                # ∫ μ(ξ) ξ dξ with μ linear on [ξ1,ξ2]
                int_μ = Δ * (μ1 + μ2) / 2
                int_μξ = Δ * (μ1 * (2ξ1 + ξ2) + μ2 * (ξ1 + 2ξ2)) / 6
                mB += int_μξ
                mA += int_μ - int_μξ
            end
            return L * mA, L * mB
        end
        msA, msB = lump_ends(μs)
        mnA, mnB = lump_ends(μn)
        jA, jB = beam_torsional_inertia ? lump_ends(μt) : (0.0, 0.0)
        if beam_axis_offsets
            pA = xa + prop.n1[1] * ŷ + prop.n1[2] * ẑ
            pB = xb + prop.n2[1] * ŷ + prop.n2[2] * ẑ
            qA = xa + prop.m1[1] * ŷ + prop.m1[2] * ẑ
            qB = xb + prop.m2[1] * ŷ + prop.m2[2] * ẑ
        else
            pA, pB, qA, qB = xa, xb, xa, xb
        end
        Jt = x̂ * x̂'
        msA > 0 && push!(lumped, LumpedMass(ga, msA, pA, jA * Jt))
        msB > 0 && push!(lumped, LumpedMass(gb, msB, pB, jB * Jt))
        mnA > 0 && push!(lumped, LumpedMass(ga, mnA, qA, Z))
        mnB > 0 && push!(lumped, LumpedMass(gb, mnB, qB, Z))
    end
    for (eid, g, cid, m, offset, J) in gfem.conm2
        p = cid == -1 ? offset : gfem.grids[g] + offset
        push!(lumped, LumpedMass(g, m, p, J))
    end
    return lumped
end

"""
    lump_by_grid(lumped) -> Vector{LumpedMass}

Combines all lumped masses attached to the same GFEM grid into a single mass m_i at x_i + δ_i with inertia J_i about x_i + δ_i, as in the methodology description.
"""
function lump_by_grid(lumped::Vector{LumpedMass})
    acc = Dict{Int,Tuple{Float64,Vector{Float64},Matrix{Float64}}}()  # m, m*p, J about origin
    for lm in lumped
        m, mp, J0 = get(acc, lm.grid, (0.0, zeros(3), zeros(3, 3)))
        acc[lm.grid] = (m + lm.m, mp + lm.m * lm.p, J0 + lm.J + transport(lm.m, lm.p))
    end
    out = LumpedMass[]
    for g in sort(collect(keys(acc)))
        m, mp, J0 = acc[g]
        m > 0 || continue
        p = mp / m
        push!(out, LumpedMass(g, m, p, J0 - transport(m, p)))
    end
    return out
end

"""
    rigid_body_mass_matrix(lumped) -> (m, cg, J_origin)

Total mass, center of gravity and inertia tensor about the basic frame origin (comparable to the Grid Point Weight Generator output).
"""
function rigid_body_mass_matrix(lumped::Vector{LumpedMass})
    m = sum(lm.m for lm in lumped)
    cg = sum(lm.m * lm.p for lm in lumped) / m
    J0 = sum(lm.J + transport(lm.m, lm.p) for lm in lumped)
    return m, cg, J0
end

# ------------------------------------------------------------------------------------------------------------
# Condensation to the beam reference axis
# ------------------------------------------------------------------------------------------------------------

"""
    straight_axis(p0, p1, nodes) -> Vector{Vector{Float64}}

Beam nodes along the straight reference axis from `p0` to `p1`. `nodes` is either the number of equally spaced nodes (end points included) or a vector of fractions in [0, 1] of the axis length.
"""
function straight_axis(p0, p1, nodes::Union{Integer,AbstractVector{<:Real}})
    ξ = nodes isa Integer ? range(0, 1, length=nodes) : nodes
    return [Float64.(p0) + s * (Float64.(p1) - Float64.(p0)) for s in ξ]
end

"""
    condense_to_beam(lumped, xref; grids, assign_by=:span, span_direction=nothing, tensor_about=:cg)

Condenses the lumped GFEM masses to the beam nodes `xref` (the reference axis) by nearest-neighbor assignment. Each GFEM grid is assigned to the beam node that is nearest in the coordinate along `span_direction` (`assign_by=:span`, i.e. the GFEM is sliced into bins centered at the beam nodes; the default direction goes from the first to the last beam node), or nearest in Euclidean distance measured from the grid position (`:grid`) or from the lumped mass position x_i + δ_i (`:mass`). If a grid is equidistant (up to round-off) from several beam nodes, its mass and inertia are split equally among them. Returns, for each beam node, the mass, the center of mass offset from the beam node, the inertia tensor about the center of mass (`tensor_about=:cg`) or about the beam node (`tensor_about=:node`), and the grid → beam node(s) assignment (a vector of nodes per grid).
"""
function condense_to_beam(lumped::Vector{LumpedMass}, xref::Vector{Vector{Float64}}; grids::Dict{Int,Vector{Float64}}, assign_by::Symbol=:span, span_direction=nothing, tensor_about::Symbol=:cg)
    n = length(xref)
    n >= 1 || throw(ArgumentError("the reference axis needs at least one node"))
    e = span_direction === nothing ? (n > 1 ? xref[end] - xref[1] : [1.0, 0.0, 0.0]) : Float64.(span_direction)
    e = e / norm(e)
    sref = [dot(e, x) for x in xref]
    m = zeros(n)
    mp = [zeros(3) for _ in 1:n]
    J0 = [zeros(3, 3) for _ in 1:n]
    # Beam nodes equidistant (up to round-off) from a grid share its lumped mass equally, instead of the floating-point argmin picking one of them
    scale = n > 1 ? maximum(norm(x - xref[1]) for x in xref) : 1.0
    tie_tol = 1e-9 * scale
    assignment = Dict{Int,Vector{Int}}()
    for lm in lumped
        nodes = get!(assignment, lm.grid) do
            d = if assign_by == :span
                s = dot(e, grids[lm.grid])
                [abs(s - sk) for sk in sref]
            elseif assign_by == :grid
                [norm(grids[lm.grid] - x) for x in xref]
            elseif assign_by == :mass
                [norm(lm.p - x) for x in xref]
            else
                throw(ArgumentError("unknown assignment criterion $assign_by"))
            end
            findall(<=(minimum(d) + tie_tol), d)
        end
        w = 1 / length(nodes)
        for k in nodes
            m[k] += w * lm.m
            mp[k] += w * lm.m * lm.p
            J0[k] += w * (lm.J + transport(lm.m, lm.p))
        end
    end
    offsets = [m[k] > 0 ? mp[k] / m[k] - xref[k] : zeros(3) for k in 1:n]
    tensors = [begin
        c = tensor_about == :cg ? xref[k] + offsets[k] : xref[k]
        J0[k] - transport(m[k], c)
    end for k in 1:n]
    return m, offsets, tensors, assignment
end

# ------------------------------------------------------------------------------------------------------------
# Output
# ------------------------------------------------------------------------------------------------------------

function write_beam_mass_model(outdir, tag, m, offsets, tensors)
    mkpath(outdir)
    open(joinpath(outdir, "beam_masses_$(tag).txt"), "w") do io
        for k in eachindex(m)
            @printf(io, "%2d %.14e\n", k, m[k])
        end
    end
    open(joinpath(outdir, "beam_offsets_$(tag).txt"), "w") do io
        for k in eachindex(m)
            @printf(io, "%2d %.14e %.14e %.14e\n", k, offsets[k]...)
        end
    end
    open(joinpath(outdir, "beam_tensors_$(tag).txt"), "w") do io
        for k in eachindex(m)
            J = tensors[k]
            @printf(io, "%2d %.14e %.14e %.14e %.14e %.14e %.14e\n", k, J[1, 1], J[2, 2], J[3, 3], J[1, 2], J[1, 3], J[2, 3])
        end
    end
end

"""
    write_beam_mass_csv(path, xref, m, offsets, tensors)

Writes the beam mass model to a CSV file with one row per beam node: node number, reference axis coordinates (x, y, z), mass, center of mass offsets (dx, dy, dz) relative to the node, and the inertia tensor components about the center of mass (Jxx, Jyy, Jzz, Jxy, Jxz, Jyz).
"""
function write_beam_mass_csv(path::AbstractString, xref, m, offsets, tensors)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "node,x,y,z,m,dx,dy,dz,Jxx,Jyy,Jzz,Jxy,Jxz,Jyz")
        for k in eachindex(m)
            J = tensors[k]
            @printf(io, "%d,%.14e,%.14e,%.14e,%.14e,%.14e,%.14e,%.14e,%.14e,%.14e,%.14e,%.14e,%.14e,%.14e\n", k, xref[k]..., m[k], offsets[k]..., J[1, 1], J[2, 2], J[3, 3], J[1, 2], J[1, 3], J[2, 3])
        end
    end
    return path
end

"""
    gpwg_from_f06(path) -> (m, cg, J_origin)

Reads the rigid-body mass matrix about the origin from the Grid Point Weight Generator output in a .f06 file.
"""
function gpwg_from_f06(path::AbstractString)
    lines = readlines(path)
    i = findfirst(l -> occursin("O U T P U T   F R O M   G R I D   P O I N T   W E I G H T", l), lines)
    i === nothing && return nothing
    MO = reduce(vcat, [parse.(Float64, split(replace(lines[i+2+k], "*" => "")))' for k in 1:6])
    m = MO[1, 1]
    # MO[1:3,4:6] = -m tilde(cg)
    cg = [MO[2, 6], MO[3, 4], MO[1, 5]] / m
    return m, cg, MO[4:6, 4:6]
end

# ------------------------------------------------------------------------------------------------------------
# Driver
# ------------------------------------------------------------------------------------------------------------

"""
    beam_mass_model(bdf, xref; shell_lumping=:rowsum, solid_lumping=:rowsum, beam_axis_offsets=true, beam_torsional_inertia=true, assign_by=:span, span_direction=nothing, tensor_about=:cg, f06=nothing, outdir=nothing, tag="beam", verbose=true)

Equivalent beam mass model of the GFEM in the Nastran input file `bdf`, condensed to the beam nodes `xref` of the user-specified reference axis (a vector of 3-vectors, see `straight_axis`; for simple prismatic models it is the elastic axis, coincident with the shear center). The reference axis is expressed in the basic coordinate system of the GFEM, in which the offsets and inertia tensors are reported. If the Nastran output `f06` is given, the lumped GFEM is validated against its Grid Point Weight Generator output. If `outdir` is given, the mass, offset and tensor files are written to it, with names ending in `tag`. The CBEAM and CBAR elements with properties in `exclude_beam_pids` and the CONM2 elements with ids in `exclude_conm2` are left out of the model. The remaining keywords are described in `lump_gfem` and `condense_to_beam`.

Returns a named tuple with the masses `m`, offsets `offsets` (center of mass relative to the beam node), inertia tensors `tensors`, grid → node(s) `assignment`, `xref`, `gfem`, `lumped` masses and the rigid-body properties `total_mass`, `cg` and `J0` (about the origin) of the lumped GFEM.
"""
function beam_mass_model(bdf::AbstractString, xref::AbstractVector; shell_lumping::Symbol=:rowsum, solid_lumping::Symbol=:rowsum, beam_axis_offsets::Bool=true, beam_torsional_inertia::Bool=true, assign_by::Symbol=:span, span_direction=nothing, tensor_about::Symbol=:cg, f06=nothing, outdir=nothing, tag::AbstractString="beam", verbose::Bool=true, exclude_beam_pids=(), exclude_conm2=())
    xref = [Float64.(x) for x in xref]
    gfem = read_gfem(bdf)
    # Elements left out of the mass model (e.g. a tip rod and its tip masses): CBEAM and CBAR elements of the properties `exclude_beam_pids` and CONM2 elements with the ids `exclude_conm2`
    filter!(b -> !(b[2] in exclude_beam_pids), gfem.beams)
    filter!(b -> !(b[2] in exclude_beam_pids), gfem.bars)
    filter!(c -> !(c[1] in exclude_conm2), gfem.conm2)
    lumped = lump_by_grid(lump_gfem(gfem; shell_lumping, solid_lumping, beam_axis_offsets, beam_torsional_inertia))
    total_mass, cg, J0 = rigid_body_mass_matrix(lumped)
    if f06 !== nothing
        gpwg = gpwg_from_f06(f06)
        if gpwg === nothing
            @warn "no Grid Point Weight Generator output found in $f06"
        elseif verbose
            @printf("GFEM lumped mass %.7e kg (GPWG %.7e), CG error %.2e m, inertia error %.2e (relative)\n", total_mass, gpwg[1], norm(cg - gpwg[2]), norm(J0 - gpwg[3]) / norm(gpwg[3]))
        end
    end
    m, offsets, tensors, assignment = condense_to_beam(lumped, xref; grids=gfem.grids, assign_by, span_direction, tensor_about)
    outdir === nothing || write_beam_mass_model(outdir, tag, m, offsets, tensors)
    return (; m, offsets, tensors, assignment, xref, gfem, lumped, total_mass, cg, J0)
end