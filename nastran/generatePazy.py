"""Generates the Nastran decks of the Pazy wing (Technion FEM, shell, CBEAM and CONM2 elements) from nastran/Pazy_Technion/S0/S0.bdf, which only has bulk data (no SPC, loads or analysis).

Usage:  python generatePazy.py [--bdf model.bdf] [--analyses static modal segments] [--nd 10] --chord-fraction 0.44 (required for static and segments) [--tag RA0p44] [--tip rigid|flexible] [--rigid-pids 13 24] --sweep 10 (required, 0 for the straight wing) [--stations y1 y2 ...] [--intervals 1] [--no-rod]

Writes next to the original model (with --tag, the static and segments files are named <name>_<tag>_<kind>.* instead):
  <name>_static.bdf    six static subcases (SOL 101), each with a unit force or moment at the load node (with --tag: <name>_<tag>_static.bdf)
  <name>_modal.bdf     normal modes (SOL 103)
  <name>_modal.json    length of the beam (--analyses modal only)
  meta.json            load node, length of the beam and other data used to post-process the static response (with --tag: <name>_<tag>_static.meta.json)
  <name>_segments.bdf  one deck with six static subcases for each window of --intervals consecutive intervals between the stations --stations: the section at its upstream station is clamped (SPC k) and the six unit loads are applied at the grid of its downstream station, which is the clamped-segment condition. The same file is overwritten for each --intervals (like for any other change of the model, e.g. --chord-fraction or --stations, the outputs of the previous Nastran run are removed and it has to be run again) (see `pazy_segments` and extractSegments.py). The whole wing is kept
  <name>_segments.json grids and positions of the stations (--analyses segments only)

Both decks are clamped at the root section of the wing (all grids at the minimum y, SPC1 in all DOF; the grids that are dependent in a rigid element, like the root grids of the RBE2 of the S10 model, are clamped by their independent grid). The beam of the stiffness model is the aluminum spar (60 x 2.25 mm, CQUAD4/CTRIA3 of property --spar-pid): its axis is at the fraction --chord-fraction (between 0 and 1) of the chord of the wing, measured from the minimum x of the shell grids, on the plane z = 0, from the root section to the end of the spar. The tip rod (CBEAM of property --rod-pid, with the tip mass at its center) is tied in the model to the wing: an RBE2 (which the model must have) from its center grid to all grids of the end section of the spar (spar and skin) makes the tip section rigid (or, as in the S10 model, one RBE2 from each of several grids of the rod, which has to be rigid, to the grids of the section next to it). In the static deck the load node is on the beam axis, at the end of the spar, and is one more dependent grid of this RBE2 (of the one whose independent grid is the closest to the load node). With --no-rod (a model without the tip rod, e.g. S0_norod.bdf) there is no rod: the load node is the independent grid of an RBE2 (all DOF) to all grids of the end section of the spar, which keeps that section rigid as the rod does, and the modal deck has no tie at the tip.

The load node uses the beam frame as its displacement and load coordinate system (CORD2R 1), so the response extracted by extractResponse.py is already in the beam frame: x1 = global y (along the span, from the root to the tip), x2 = -global x (chord, towards the leading edge, if the model has x aft) and x3 = global z (thickness), which is right-handed.

With --sweep (angle in degrees, positive towards the trailing edge, for models like S10 and S20 whose sections, in planes of constant y, are shifted along x by y tan(sweep)) the beam axis is inclined by that angle, at a constant fraction of the chord along the span (the x of the sections and of the chord extent are taken without the shift, x - y tan(sweep)), and the beam frame is rotated by the angle about z so that x1 is along the axis: x1 = (sin(sweep), cos(sweep), 0), x2 = (-cos(sweep), sin(sweep), 0) and x3 = global z. The length of the beam is the length of the axis, (y of the end of the spar - y of the root)/cos(sweep). The segments deck is not available for swept wings (the generalized displacements of the sections can be taken from the static job, with extractStations.py).
"""
import argparse
import json
from pathlib import Path

import numpy as np
from pyNastran.bdf.bdf import BDF
from pyNastran.bdf.case_control_deck import CaseControlDeck

from common import NASTRAN_DIR, write_bdf_if_changed, write_meta

DEFAULT_BDF = NASTRAN_DIR / "Pazy_Technion" / "S0" / "S0.bdf"
SPAR_PID = 100003
ROD_PID = 13
BEAM_FRAME = 1
RIGID_MID = 100005  # material of the rigid tip rod (--tip rigid)
RIGID_MODULUS = 1e15  # Young's modulus [Pa] of the rigid tip rod
TIP_SUFFIXES = {"rigid": "_rigidtip", "flexible": "_flextip"}  # added to the names of the decks and jobs of --tip
TOL = 1e-6
STATION_TOL = 1e-5
UNIT_LOADS = [
    ("FORCE", [1., 0., 0.]), ("FORCE", [0., 1., 0.]), ("FORCE", [0., 0., 1.]),
    ("MOMENT", [1., 0., 0.]), ("MOMENT", [0., 1., 0.]), ("MOMENT", [0., 0., 1.]),
]


def add_unit_load(model, sid, kind, grid, direction):
    """Unit force (N) or moment (N*m) at `grid` in the beam frame."""
    if kind == "FORCE":
        model.add_force(sid, grid, 1.0, direction, cid=BEAM_FRAME)
    else:
        model.add_moment(sid, grid, 1.0, direction, cid=BEAM_FRAME)


def read_model(bdf, rod_pid=None, rigid_pids=()):
    """Model of the bulk data `bdf`. With `rigid_pids`, the CBEAM properties of the tip rod (the properties that the model lacks are skipped) get a new material, rigid (Young's modulus RIGID_MODULUS) and with the Poisson's ratio and density of the material of the property `rod_pid`, i.e. a rigid nylon, as the one of the tip rod and breaking beam of the S10 and S20 models."""
    model = BDF(debug=False)
    model.read_bdf(str(bdf), xref=False, punch=True)
    pids = [pid for pid in rigid_pids if pid in model.properties]
    if pids:
        if RIGID_MID in model.materials:
            raise SystemExit(f"Material {RIGID_MID} already exists")
        mat = model.materials[model.properties[rod_pid].mid]
        model.add_mat1(RIGID_MID, RIGID_MODULUS, None, mat.nu, rho=mat.rho)
        for pid in pids:
            model.properties[pid].mid = RIGID_MID
    return model


def tip_rod_center(model, rod_pid):
    """Grid at the center of the tip rod (closest to the mean position of its grids) and the grids of the rod. Without rod (`rod_pid` None): (None, [])."""
    if rod_pid is None:
        return None, []
    rod = sorted({n for e in model.elements.values() if e.type == "CBEAM" and e.pid == rod_pid for n in e.node_ids})
    if not rod:
        raise SystemExit(f"No CBEAM elements with property {rod_pid} in the model")
    rod_xyz = np.array([model.nodes[i].xyz for i in rod])
    center = rod[int(np.argmin(np.linalg.norm(rod_xyz - rod_xyz.mean(axis=0), axis=1)))]
    return center, rod


def beam_geometry(model, spar_pid):
    """Root section grids (all grids at the minimum y, except those that are dependent in a rigid element, which are clamped by their independent grid), grids of the end section of the spar (spar and skin, without the tip rod): (root nodes, tip nodes, y of the root, y of the spar end). The beam axis is set by `beam_axis_x`."""
    ids = np.array(list(model.nodes))
    xyz = np.array([model.nodes[i].xyz for i in ids])
    y_root = xyz[:, 1].min()
    dependent = {n for r in model.rigid_elements.values() for n in r.Gmi_node_ids}
    root_nodes = [n for n in ids[np.abs(xyz[:, 1] - y_root) < TOL].tolist() if n not in dependent]

    spar = sorted({n for e in model.elements.values() if e.type in ("CQUAD4", "CTRIA3") and e.pid == spar_pid for n in e.node_ids})
    if not spar:
        raise SystemExit(f"No shell elements with property {spar_pid} in the model")
    spar_xyz = np.array([model.nodes[i].xyz for i in spar])
    y_end = spar_xyz[:, 1].max()
    tip_nodes = ids[np.abs(xyz[:, 1] - y_end) < TOL].tolist()
    return root_nodes, tip_nodes, y_root, y_end


def beam_axis_x(model, chord_fraction, sweep):
    """x of the beam axis at y = 0, at the fraction `chord_fraction` (between 0 and 1) of the chord, measured from the minimum x of the grids of all shell elements. For a swept wing (`sweep` in degrees, with the sections in planes of constant y shifted along x by y tan(sweep)) the chordwise positions are those without the shift, x - y tan(sweep)."""
    if not 0 <= chord_fraction <= 1:
        raise ValueError(f"chord_fraction must be between 0 and 1, got {chord_fraction}")
    tan_sweep = np.tan(np.radians(sweep))
    shell_x = np.array([model.nodes[n].xyz[0] - tan_sweep * model.nodes[n].xyz[1] for e in model.elements.values() if e.type in ("CQUAD4", "CTRIA3") for n in e.node_ids])
    return shell_x.min() + chord_fraction * (shell_x.max() - shell_x.min())


def beam_length(y_root, y_end, sweep):
    """Length of the beam axis, from the root section to the end of the spar (`sweep` in degrees)."""
    return float((y_end - y_root) / np.cos(np.radians(sweep)))


def tip_rod_rbe2(model, rod_nodes, tip_nodes, point=None):
    """RBE2 of the model that ties the tip rod to the grids of the tip section (all DOF), to which the load node is added. The tip section is tied to the rod by one RBE2 from the center of the rod to all of its grids (Technion model) or by one RBE2 from each of several grids of the rod (which is rigid) to the grids of the section next to it (S10 model). Returns the RBE2 whose independent grid is the closest to `point` (default: the center of the tip section)."""
    ties = [r for r in model.rigid_elements.values() if r.type == "RBE2" and set(r.independent_nodes) & set(rod_nodes)]
    tied = {n for r in ties for n in r.Gmi_node_ids}
    if not ties or not tied <= set(tip_nodes) or any(len(r.independent_nodes) != 1 or str(r.cm) != "123456" for r in ties):
        raise SystemExit(f"The model must have RBE2 (all DOF) from the grids of the tip rod only to the grids of the end section of the spar ({len(tip_nodes)} grids)")
    if point is None:
        point = np.mean([model.nodes[n].xyz for n in tip_nodes], axis=0)
    return min(ties, key=lambda r: np.linalg.norm(np.asarray(model.nodes[r.independent_nodes[0]].xyz) - point))


def tip_load_rbe2(model, load_node, tip_nodes):
    """RBE2 (all DOF) of a model without tip rod that ties the load node, its independent grid, to the grids of the end section of the spar, so that the tip section is rigid as it is with the rod."""
    eid = max(max(model.rigid_elements, default=0), max(model.elements)) + 1
    model.add_rbe2(eid, load_node, "123456", tip_nodes)


def spar_section_grids(model, spar_pid, y, tol=STATION_TOL):
    """Grids of the spar (shell elements of property `spar_pid`) at the spanwise position `y`, within `tol`."""
    spar = sorted({n for e in model.elements.values() if e.type in ("CQUAD4", "CTRIA3") and e.pid == spar_pid for n in e.node_ids})
    return [n for n in spar if abs(model.nodes[n].xyz[1] - y) < tol]


def station_kind(y, y_root, y_end, tol=STATION_TOL):
    """Station at the clamped root, at the end of the spar (the load node) or in between."""
    if abs(y - y_root) < tol:
        return "root"
    if abs(y - y_end) < tol:
        return "tip"
    return "node"


def rotation_weights(model, section, x_rot, tol=1e-9):
    """Spar grids of a section and their weights to interpolate linearly, at the chordwise position `x_rot`, a quantity at the grids: the two grids that enclose it (a single one if it is at a grid)."""
    xs = np.array([model.nodes[n].xyz[0] for n in section])
    order = np.argsort(xs)
    section, xs = [section[i] for i in order], xs[order]
    if not xs[0] - tol <= x_rot <= xs[-1] + tol:
        raise SystemExit(f"x = {x_rot} is outside the chord of the spar grids ({xs[0]}, {xs[-1]})")
    i = int(np.clip(np.searchsorted(xs, x_rot) - 1, 0, len(xs) - 2))
    wb = (x_rot - xs[i]) / (xs[i + 1] - xs[i])
    if wb < tol:
        return [section[i]], [1.0]
    if wb > 1 - tol:
        return [section[i + 1]], [1.0]
    return [section[i], section[i + 1]], [round(1 - wb, 8), round(wb, 8)]


def add_station_grids(model, stations, spar_pid, x_axis, y_root, y_end, load_node, tol=STATION_TOL):
    """Grid on the beam axis at each intermediate station, with the beam frame as displacement coordinate system, tied by an RBE3 to the spar grids at the same y. The RBE3 does not add stiffness. The translations, and the rotations about the span and normal axes, are the rigid body motion fitted to the translations of all grids. The rotation about the chord axis, which the translations of grids on the mid-plane cannot show, is the rotation of the grids about x (the normal of the plate) linearly interpolated at the chordwise position of the axis between the two grids that enclose it, instead of the mean of the rotations of all grids, which would be the value at their centroid: the normal of the plate does not rotate rigidly with the section, but varies along the chord with the twist rate. Returns the list of stations as dicts with the position `y` and the `grid` of the station (the load node at the end of the spar, None at the root, which is clamped)."""
    dependent = {n for r in model.rigid_elements.values() for n in r.Gmi_node_ids}
    out = []
    for y in stations:
        kind = station_kind(y, y_root, y_end, tol)
        grid = {"root": None, "tip": load_node}.get(kind)
        if kind == "node":
            section = spar_section_grids(model, spar_pid, y, tol)
            if not section:
                raise SystemExit(f"No spar grids at y = {y} (within {tol}): the stations must be at rows of grids of the spar")
            if dependent & set(section):
                raise SystemExit(f"Spar grids at y = {y} are dependent in other rigid elements: {sorted(dependent & set(section))}")
            enclosing, weights = rotation_weights(model, section, x_axis)
            grid = max(model.nodes) + 1
            model.add_grid(grid, [x_axis, y, 0.], cd=BEAM_FRAME)
            eid = max(max(model.rigid_elements, default=0), max(model.elements)) + 1
            model.add_rbe3(eid, grid, "123456", [1.0] + weights, ["123"] + ["4"] * len(weights), [section] + [[n] for n in enclosing])
        out.append(dict(y=float(y), grid=None if grid is None else int(grid), kind=kind))
    return out


def _finish(model, case_lines, filename, sol):
    model.punch = False  # the original deck has bulk data only; the new one is complete
    model.sol = sol
    model.case_control_deck = CaseControlDeck(case_lines)
    model.add_param("POST", [-1])  # generates .op2
    write_bdf_if_changed(model, filename, size=16)


def prepare_wing(bdf, spar_pid, rod_pid, chord_fraction, sweep, rigid_pids=()):
    """Wing model with the beam frame, the load node at the end of the beam axis and the tip rod tied to the tip section (no constraints or loads yet). Returns the model, the load node, the root grids, the rod grids and the geometry: (model, load_node, root_nodes, tip_nodes, rod_nodes, rod_center, y_root, y_end, x_axis), with x_axis the position of the axis at y = 0 (see `beam_axis_x`)."""
    model = read_model(bdf, rod_pid, rigid_pids)
    root_nodes, tip_nodes, y_root, y_end = beam_geometry(model, spar_pid)
    x_axis = beam_axis_x(model, chord_fraction, sweep)
    rod_center, rod_nodes = tip_rod_center(model, rod_pid)
    if set(tip_nodes) & set(rod_nodes):
        raise SystemExit("The tip rod is in the end section of the spar")

    # Load node on the beam axis at the end of the spar, coincident with the last node of the beam mass model, with the beam frame as coordinate system (the basic frame has x = chord, y = span, z = thickness): x1 along the axis, which is inclined by `sweep` from y, and x3 = z
    s, c = np.sin(np.radians(sweep)), np.cos(np.radians(sweep))
    x_end = x_axis + y_end * np.tan(np.radians(sweep))
    load_node = max(model.nodes) + 1
    model.add_cord2r(BEAM_FRAME, [x_end, y_end, 0.], [x_end, y_end, 1.], [x_end + s, y_end + c, 0.], rid=0)
    model.add_grid(load_node, [x_end, y_end, 0.], cd=BEAM_FRAME)

    # The load node is rigidly tied to the tip rod together with the grids of the tip section (one more dependent grid of an RBE2 of the model), or without rod to the grids of the tip section by an RBE2 of its own
    if rod_nodes:
        tip_rod_rbe2(model, rod_nodes, tip_nodes, [x_end, y_end, 0.]).Gmi.append(load_node)
    else:
        tip_load_rbe2(model, load_node, tip_nodes)
    return model, load_node, root_nodes, tip_nodes, rod_nodes, rod_center, y_root, y_end, x_axis


def pazy_static(bdf, filename, chord_fraction, sweep, spar_pid=SPAR_PID, rod_pid=ROD_PID, rigid_pids=()):
    """Six static subcases (SOL 101), each with one unit force or moment (N and N*m) at the load node, with the root section clamped (SPC 1). `chord_fraction` places the beam axis (see `beam_axis_x`) and `sweep` is the sweep angle in degrees."""
    model, load_node, root_nodes, tip_nodes, rod_nodes, rod_center, y_root, y_end, x_axis = prepare_wing(bdf, spar_pid, rod_pid, chord_fraction, sweep, rigid_pids)

    model.add_spc1(1, "123456", root_nodes)

    case_lines = ["ECHO = NONE"]
    for sid, (kind, direction) in enumerate(UNIT_LOADS, start=1):
        add_unit_load(model, sid, kind, load_node, direction)
        case_lines += [
            f"SUBCASE {sid}",
            f"  LOAD = {sid}",
            "  SPC = 1",
            "  DISPLACEMENT(PLOT) = ALL",
        ]
    _finish(model, case_lines, filename, sol=101)

    meta = dict(load_node=load_node, L=beam_length(y_root, y_end, sweep), y_root=float(y_root), y_end=float(y_end), x_axis=float(x_axis), n_root_nodes=len(root_nodes), n_tip_nodes=len(tip_nodes), rod_center=None if rod_center is None else int(rod_center), spar_pid=spar_pid, rod_pid=rod_pid)
    meta["sweep"] = float(sweep)
    return meta


def pazy_segments(bdf, filename, stations, chord_fraction, sweep, spar_pid=SPAR_PID, rod_pid=ROD_PID, station_tol=STATION_TOL, intervals=1, rigid_pids=()):
    """One deck with six static subcases (SOL 101) for each window w of `intervals` consecutive intervals between stations (w = 1 ... n - intervals, n stations): the section of the wing at station w is clamped (SPC w, all the grids of the model at that y, as the root section) and the six unit loads (N and N*m) are applied at the grid of the station w + intervals (see `add_station_grids`; the load node at the end of the spar). The rest of the wing is kept, free and unloaded: upstream of the clamp it carries no load, and downstream it stiffens the end of the window as in the wing. Subcases 6(w-1)+1 ... 6w are those of the window w, with load cases [F1, F2, F3, M1, M2, M3]. The loads are applied where the displacements are taken (the RBE3 is symmetric in the work of the generalized forces), so the displacement-force matrix of the window is reciprocal. With more than one interval the end effects of the clamp and of the load (warping) are spread over a longer beam. `stations` must go from the root to the end of the spar."""
    model, load_node, root_nodes, tip_nodes, rod_nodes, rod_center, y_root, y_end, x_axis = prepare_wing(bdf, spar_pid, rod_pid, chord_fraction, sweep, rigid_pids=rigid_pids)
    n = len(stations)
    if intervals < 1 or n < intervals + 1:
        raise SystemExit(f"{intervals} intervals per window need at least {intervals + 1} stations, got {n}")

    # Grids of the original model to clamp at each station, found before the grids of the stations (at the same y) are added
    ids = [n for n in model.nodes if n not in rod_nodes and n != load_node]
    ys = np.array([model.nodes[n].xyz[1] for n in ids])
    clamped = []
    for k, y in enumerate(stations[:n - intervals], start=1):
        if station_kind(y, y_root, y_end, station_tol) == "tip":
            raise SystemExit("The stations must be increasing, from the root to the end of the spar")
        if k == 1 and station_kind(y, y_root, y_end, station_tol) != "root":
            raise SystemExit(f"The first station (y = {y}) must be the root section (y = {y_root})")
        section = root_nodes if k == 1 else [n for n, yn in zip(ids, ys) if abs(yn - y) < station_tol]
        if not section:
            raise SystemExit(f"No grids at y = {y} (within {station_tol})")

        # Grids that are dependent in a rigid element cannot be constrained, but they are clamped with their independent grid, which has to be in the section
        grids = set(section)
        for r in model.rigid_elements.values():
            if grids & set(r.Gmi_node_ids) and not set(r.independent_nodes) <= grids:
                raise SystemExit(f"Rigid element {r.eid} ties grids at y = {y} to grids out of that section")
        grids = sorted(grids - {n for r in model.rigid_elements.values() for n in r.Gmi_node_ids})
        model.add_spc1(k, "123456", grids)
        clamped.append(len(grids))

    station_grids = add_station_grids(model, stations, spar_pid, x_axis, y_root, y_end, load_node, station_tol)
    if station_grids[-1]["kind"] != "tip":
        raise SystemExit(f"The last station (y = {stations[-1]}) must be the end of the spar (y = {y_end})")

    # Displacements only at the grids of the stations (they are the loaded ones), to keep the .op2 small
    case_lines = ["ECHO = NONE", "SET 1 = " + ", ".join(str(s["grid"]) for s in station_grids if s["grid"] is not None)]
    for k in range(1, n - intervals + 1):
        for j, (kind, direction) in enumerate(UNIT_LOADS, start=1):
            sid = 6 * (k - 1) + j
            add_unit_load(model, sid, kind, station_grids[k - 1 + intervals]["grid"], direction)
            case_lines += [
                f"SUBCASE {sid}",
                f"  LOAD = {sid}",
                f"  SPC = {k}",
                "  DISPLACEMENT(PLOT) = 1",
            ]
    _finish(model, case_lines, filename, sol=101)

    return dict(load_node=load_node, L=float(y_end - y_root), y_root=float(y_root), y_end=float(y_end), x_axis=float(x_axis), n_clamped_grids=clamped, rod_center=None if rod_center is None else int(rod_center), spar_pid=spar_pid, rod_pid=rod_pid, stations=station_grids, intervals=intervals)


def pazy_modal(bdf, filename, sweep, nd=10, spar_pid=SPAR_PID, rod_pid=ROD_PID, rigid_pids=()):
    """Normal modes (SOL 103) of the wing clamped at the root section, with the tip section tied to the tip rod. nd: number of modes requested. `sweep` (degrees) only sets the length of the beam in the returned data (the length of the axis, see `beam_length`)."""
    model = read_model(bdf, rod_pid, rigid_pids)
    root_nodes, tip_nodes, y_root, y_end = beam_geometry(model, spar_pid)
    _, rod_nodes = tip_rod_center(model, rod_pid)
    if rod_nodes:
        tip_rod_rbe2(model, rod_nodes, tip_nodes)
    model.add_spc1(1, "123456", root_nodes)
    model.add_eigrl(10, nd=nd, norm="MASS")
    case_lines = [
        "ECHO = NONE",
        "SUBCASE 1",
        "  SPC = 1",
        "  METHOD = 10",
        "  DISPLACEMENT(PLOT) = ALL",
    ]
    _finish(model, case_lines, filename, sol=103)
    return dict(L=beam_length(y_root, y_end, sweep), y_root=float(y_root), y_end=float(y_end), sweep=float(sweep))


def unit_interval(value):
    x = float(value)
    if not 0 <= x <= 1:
        raise argparse.ArgumentTypeError(f"must be between 0 and 1, got {value}")
    return x


def main():
    parser = argparse.ArgumentParser(description="Generate the static, modal and segments BDF files of the Pazy wing next to the original model")
    parser.add_argument("--bdf", type=Path, default=DEFAULT_BDF, help="Pazy wing model with bulk data only")
    parser.add_argument("--analyses", nargs="+", choices=["static", "modal", "segments"], default=["static", "modal"])
    parser.add_argument("--nd", type=int, default=10, help="number of modes requested")
    parser.add_argument("--spar-pid", type=int, default=SPAR_PID, help="property of the shell elements of the spar")
    parser.add_argument("--chord-fraction", type=unit_interval, default=None, help="chordwise position of the beam axis and load node as a fraction (between 0 and 1) of the chord from the minimum x of the shell grids (required for the static and segments analyses)")
    parser.add_argument("--sweep", type=float, required=True, help="sweep angle [deg] of the wing (0 for a straight wing), whose sections (planes of constant y) are shifted along x by y tan(sweep), e.g. the S10 and S20 models: inclines the beam axis and the beam frame (the segments analysis is only available with 0)")
    parser.add_argument("--tag", help="label added to the names of the static and segments decks, <name>_<tag>_static.bdf (which then has its own <name>_<tag>_static.meta.json), so that decks with different beam axes (--chord-fraction) or stations are kept side by side. The modal deck does not depend on them and is not tagged")
    parser.add_argument("--rod-pid", type=int, default=ROD_PID, help="property of the CBEAM elements of the tip rod")
    parser.add_argument("--tip", choices=sorted(TIP_SUFFIXES), help="flexible or rigid tip rod, for the models with the rod: the names of the decks and jobs get the suffix _flextip or _rigidtip (e.g. S0_RA0p43_static_rigidtip.bdf, S0_modal_flextip.bdf, so that the results of both are kept), and with `rigid` the material of the rod is a rigid one (Young's modulus 1e15 Pa, with the Poisson's ratio and density of the original), so that the tip section of the wing is rigid for the extraction of the stiffness. With `flexible` the material is the one of the model")
    parser.add_argument("--rigid-pids", nargs="+", type=int, help="properties that get the rigid material with --tip rigid (default: --rod-pid; the swept models S10 and S20 also have the breaking beam, 24, in LE and TE). Those that the model lacks are skipped, and the Poisson's ratio and density of the rigid material are those of the material of --rod-pid")
    parser.add_argument("--no-rod", action="store_true", help="the model has no tip rod (e.g. S0_norod.bdf): the load node is tied to the tip section by an RBE2 of its own (static and segments decks), and the modal deck has no tie at the tip")
    parser.add_argument("--stations", nargs="+", type=float, help="spanwise positions (y) of the nodes of the discretized beam, from the root to the end of the spar, for the segments deck")
    parser.add_argument("--tol", type=float, default=STATION_TOL, help="tolerance in y to match the stations with the rows of grids of the spar")
    parser.add_argument("--intervals", type=int, default=1, help="number of intervals between stations of each clamped window, segments deck (1, 2, 4, ...)")
    args = parser.parse_args()
    rod_pid = None if args.no_rod else args.rod_pid
    if args.sweep and "segments" in args.analyses:
        parser.error("the segments analysis is not available with --sweep")
    if args.chord_fraction is None and {"static", "segments"} & set(args.analyses):
        parser.error("--chord-fraction is required for the static and segments analyses")

    bdf = args.bdf.resolve()
    directory = bdf.parent
    name = bdf.stem

    tip_suffix = TIP_SUFFIXES[args.tip] if args.tip else ""
    rigid_pids = tuple(args.rigid_pids if args.rigid_pids is not None else [args.rod_pid]) if args.tip == "rigid" and rod_pid is not None else ()

    def job_name(kind):
        return (f"{name}_{args.tag}_{kind}" if args.tag else f"{name}_{kind}") + tip_suffix

    if "modal" in args.analyses:
        modal_file = directory / f"{name}_modal{tip_suffix}.bdf"  # not tagged
        meta = pazy_modal(bdf, modal_file, sweep=args.sweep, nd=args.nd, spar_pid=args.spar_pid, rod_pid=rod_pid, rigid_pids=rigid_pids)
        # Written after the BDF, which removes the stale outputs of the job (same stem) when it changes
        modal_file.with_suffix(".json").write_text(json.dumps(dict(model=bdf.name, **meta), indent=2))
    if "static" in args.analyses:
        static_name = job_name("static")
        meta = pazy_static(bdf, directory / f"{static_name}.bdf", spar_pid=args.spar_pid, rod_pid=rod_pid, chord_fraction=args.chord_fraction, sweep=args.sweep, rigid_pids=rigid_pids)
        # Written after the BDF, which removes the stale outputs of the job (same stem) when it changes
        write_meta(directory, job=static_name if args.tag or tip_suffix else None, model=bdf.name, **meta)
        print("Load node:", meta["load_node"])
        print("Length of the beam [m]:", meta["L"])
    if "segments" in args.analyses:
        if not args.stations:
            parser.error("--stations is required for the segments analysis")
        segments_file = directory / f"{job_name('segments')}.bdf"
        meta = pazy_segments(bdf, segments_file, args.stations, spar_pid=args.spar_pid, rod_pid=rod_pid, chord_fraction=args.chord_fraction, sweep=args.sweep, station_tol=args.tol, intervals=args.intervals, rigid_pids=rigid_pids)
        segments_file.with_suffix(".json").write_text(json.dumps(dict(model=bdf.name, **meta), indent=2))
        windows = len(args.stations) - args.intervals
        print("Windows of", args.intervals, "intervals:", windows, "with", 6 * windows, "subcases")
    print("BDF files in", directory)


if __name__ == "__main__":
    main()
