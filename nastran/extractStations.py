"""Extracts the generalized displacements of the nodes of the discretized beam (stations) of the Pazy wing from the Nastran .op2 of the six static load cases at the end of the beam.

Usage:  python extractStations.py "nastran/Pazy_Technion/S0" --job S0_static --stations y1 y2 ... [--twist-center axis|auto|<chord fraction>]

The response is the rigid body motion of the section, fitted to the displacements of all spar grids (property --spar-pid) at the same y (within --tol), referred to the point of the beam axis (see `section_motion`). The rotation about the chord axis is that of the plate at the twist center (--twist-center, see `twist_center`), so that the response does not depend on where the beam axis is. Any static job can be used, because the .op2 has the displacements of all grids; the axis is in meta.json (written by generatePazy.py --analyses static) and --stations are the spanwise positions of the nodes. For a swept wing (sweep in meta.json, generatePazy.py --sweep), the sections are still those of constant y, the point of the axis at the position y is at x_axis + y tan(sweep) and the response is in the beam frame rotated by the sweep angle (as that of the load node); the station is at the distance (y - y_root)/cos(sweep) from the root along the axis.

Writes <job>.station_response_section.csv (<job>.station_response_section_<tag>.csv with --tag): the stations stacked, (6 n x 6) matrix where rows 6(k-1)+1:6k are [Tx, Ty, Tz, Rx, Ry, Rz] of station k (m and rad, in the beam frame: x1 = y of the Nastran model along the span, x2 = -x, x3 = z; with sweep, x1 = (sin(sweep), cos(sweep), 0), x2 = (-cos(sweep), sin(sweep), 0), x3 = z) and columns are the load cases (Fx, Fy, Fz, Mx, My, Mz). The root station is clamped (zeros).
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2

from common import read_meta
from generatePazy import beam_axis_x

SPAR_PID = 100003
TOL = 1e-5
MIN_SECTION_GRIDS = 3  # minimum number of spar grids at a station to fit its rigid body motion


def meta_sweep(meta, job):
    """Sweep angle [deg] that generatePazy.py wrote in the metadata of the job."""
    if "sweep" not in meta:
        sys.exit(f"The metadata of {job} has no sweep angle: generate the deck again with generatePazy.py (--sweep)")
    return meta["sweep"]


def basic_to_beam(sweep):
    """Rotation from the basic frame of the Nastran model (x = chord, y = span, z = thickness) to the beam frame, whose rows are x1 = (sin(sweep), cos(sweep), 0) along the axis, x2 = (-cos(sweep), sin(sweep), 0) and x3 = z (`sweep` in degrees; without sweep x1 = y, x2 = -x)."""
    s, c = np.sin(np.radians(sweep)), np.cos(np.radians(sweep))
    return np.array([[s, c, 0.], [-c, s, 0.], [0., 0., 1.]])


def read_displacements(op2_file):
    """Displacements of the six load cases: (grid ids, array (6, n grids, 6))."""
    op2 = OP2(debug=False)
    op2.read_op2(str(op2_file))
    ids = op2.displacements[1].node_gridtype[:, 0]
    return ids, np.array([op2.displacements[sid].data[0] for sid in range(1, 6 + 1)], dtype=float)


def skew(r):
    return np.array([[0., -r[2], r[1]], [r[2], 0., -r[0]], [-r[1], r[0], 0.]])


def section_motion(points, motion, p0, x_rot=None):
    """Rigid body motion of the section with spar grids at `points` (m, 3) and generalized displacements `motion` (m, 6), in the basic frame, referred to the point `p0` of the beam axis. Returns the (6,) vector [u; θ].

    The translation u and the rotations about the span (y) and the normal (z) are the least squares rigid body motion of the translations of the grids. The rotation about the chord axis (x) cannot be seen in the translations of grids on the mid-plane, so it comes from the rotations of the grids about x, which are those of the normal of the plate, -∂w/∂y. The section does not rotate rigidly: the slope of the plate varies linearly along the chord with the twist rate (warping), θx(x) = θ - κ (x - x_tc), where the warping vanishes at the twist center x_tc. So the rigid body rotation of the section is the value at x_tc of the linear fit of the rotations of the grids, `x_rot` (default: the chordwise position of `p0`, which is only right if the axis is at the twist center; at any other point the result mixes the warping of the section into the rotation, and the beam model is not reciprocal)."""
    r = points - p0
    # Translation of the grid: u + θ x r = u - skew(r) θ. θx does not move grids on the line of the section
    A = np.vstack([np.hstack([np.eye(3), -skew(ri)])[:, [0, 1, 2, 4, 5]] for ri in r])
    x, *_ = np.linalg.lstsq(A, motion[:, :3].reshape(-1), rcond=None)
    x_rot = p0[0] if x_rot is None else x_rot
    theta_x = np.polyfit(points[:, 0] - x_rot, motion[:, 3], 1)[1]
    return np.array([x[0], x[1], x[2], theta_x, x[3], x[4]])


def to_beam_frame(v, sweep):
    R = basic_to_beam(sweep)
    return np.concatenate([R @ v[:3], R @ v[3:]])


def spar_grids(model, spar_pid):
    return sorted({n for e in model.elements.values() if e.type in ("CQUAD4", "CTRIA3") and e.pid == spar_pid for n in e.node_ids})


def spar_rows(model, spar, tol):
    """Spanwise positions of the rows of spar grids (grids whose y differ by less than `tol` are in the same row)."""
    ys = np.sort([model.nodes[n].xyz[1] for n in spar])
    rows = [[ys[0]]]
    for y in ys[1:]:
        if y - rows[-1][-1] > tol:
            rows.append([])
        rows[-1].append(y)
    return np.array([np.mean(r) for r in rows])


def twist_center(model, spar, rows, disp, x_axis, y_root, y_end, tol, sweep, margin=0.2):
    """x of the twist center of the spar (at y = 0, like `x_axis`, so for a swept wing it is shifted by y tan(sweep) at each section), from the reciprocity of the beam model: the compliance of the rotations under the moments at the end of the beam is the integral of the (symmetric) block C22 of the compliance, so the rotation about the normal axis under the moment about the span, θ2(M1), and the rotation about the span under the moment about the normal axis, θ1(M2), of a section must be equal. θ2 is the rotation about the chord axis, which depends linearly on the chordwise position at which it is taken (see `section_motion`), so each row of spar grids gives the position at which they are equal. Only the rows in the central part of the span are used, away from the end effects of the clamp and of the tip (`margin` is the fraction of the span discarded at each end). `rows` maps the grid ids to their index in `disp`. Returns the x of the twist center, the array with the values of the rows and the y of the rows."""
    tan_sweep = np.tan(np.radians(sweep))
    estimates, ys = [], []
    for y in spar_rows(model, spar, tol):
        if not y_root + margin * (y_end - y_root) <= y <= y_end - margin * (y_end - y_root):
            continue
        section = [n for n in spar if abs(model.nodes[n].xyz[1] - y) < tol]
        points = np.array([model.nodes[n].xyz for n in section])
        p0 = np.array([x_axis + y * tan_sweep, y, 0.])
        idx = [rows[n] for n in section]

        def mismatch(x_rot):
            theta2_M1 = to_beam_frame(section_motion(points, disp[3, idx, :], p0, x_rot), sweep)[4]
            theta1_M2 = to_beam_frame(section_motion(points, disp[4, idx, :], p0, x_rot), sweep)[3]
            return theta2_M1 - theta1_M2

        # The mismatch is linear in the position of the rotation
        xa, xb = p0[0], p0[0] + 1e-3
        da, db = mismatch(xa), mismatch(xb)
        if db != da:
            estimates.append(xa - da * (xb - xa) / (db - da) - y * tan_sweep)
            ys.append(y)
    if not estimates:
        sys.exit("The twist center cannot be found: no rows of spar grids in the central part of the span with a rotation that depends on the chordwise position")
    return float(np.median(estimates)), np.array(estimates), np.array(ys)


def section_response(ids, disp, bdf, stations, x_axis, y_root, y_end, spar_pid, tol, sweep, twist_center_arg="axis"):
    """Generalized displacements of the stations, (n stations, 6, 6). `twist_center_arg` is the chordwise position at which the rotation about the chord axis of the sections is taken (see `section_motion`): "axis" (that of the beam axis), "auto" (found from the response, see `twist_center`) or a fraction of the chord (see `beam_axis_x`)."""
    model = BDF(debug=False)
    model.read_bdf(str(bdf), xref=False, punch=True)
    spar = spar_grids(model, spar_pid)
    if any(model.nodes[n].cd != 0 for n in spar):
        sys.exit("Spar grids with a non-basic displacement coordinate system are not supported")
    rows = {int(g): i for i, g in enumerate(ids)}
    tan_sweep = np.tan(np.radians(sweep))

    # Chordwise position of the rotation of the sections, at y = 0
    x0, c = beam_axis_x(model, 0.0, sweep), beam_axis_x(model, 1.0, sweep) - beam_axis_x(model, 0.0, sweep)
    if twist_center_arg == "axis":
        x_tc = None
    elif twist_center_arg == "auto":
        x_tc, estimates, ys = twist_center(model, spar, rows, disp, x_axis, y_root, y_end, tol, sweep)
        print(f"Twist center (from {len(estimates)} rows of spar grids, y from {ys.min():.4f} to {ys.max():.4f}): x = {x_tc:.6f} m, {(x_tc - x0) / c:.4f} of the chord (rows: {(estimates.min() - x0) / c:.4f} to {(estimates.max() - x0) / c:.4f})")
    else:
        x_tc = beam_axis_x(model, float(twist_center_arg), sweep)
        print(f"Twist center: x = {x_tc:.6f} m, {float(twist_center_arg):.4f} of the chord")

    T = np.zeros((len(stations), 6, 6))
    for k, y in enumerate(stations):
        if abs(y - y_root) < tol:
            continue  # clamped
        section = [n for n in spar if abs(model.nodes[n].xyz[1] - y) < tol]
        if not section:
            sys.exit(f"No spar grids at y = {y} (within {tol}): the stations must be at rows of grids of the spar")
        if len(section) < MIN_SECTION_GRIDS:
            sys.exit(f"Only {len(section)} spar grid(s) at y = {y} (within {tol}): the rigid body motion of the section is not determined with fewer than {MIN_SECTION_GRIDS}; use a row of grids that spans the chord of the spar")
        points = np.array([model.nodes[n].xyz for n in section])
        p0 = np.array([x_axis + y * tan_sweep, y, 0.])
        x_rot = None if x_tc is None else x_tc + y * tan_sweep
        for j in range(6):
            motion = disp[j, [rows[n] for n in section], :]
            T[k, :, j] = to_beam_frame(section_motion(points, motion, p0, x_rot), sweep)
        print(f"y = {y:.6f} m: {len(section)} spar grids")
    return T


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--job", required=True, help="name of the static job, without extension (e.g. S0_static)")
    parser.add_argument("--stations", nargs="+", type=float, required=True, help="spanwise positions (y) of the nodes of the discretized beam")
    parser.add_argument("--spar-pid", type=int, default=SPAR_PID, help="property of the shell elements of the spar")
    parser.add_argument("--tag", help="label added to the name of the output, <job>.station_response_section_<tag>.csv, e.g. to keep the responses of different stations of the same job")
    parser.add_argument("--tol", type=float, default=TOL, help="tolerance in y to take the grids of the spar at a station")
    parser.add_argument("--twist-center", default="axis", help="chordwise position at which the rotation about the chord axis of the sections is taken: 'axis' (the beam axis, which gives a response that depends on the axis), 'auto' (the twist center, found from the reciprocity of the response, so that it does not depend on the axis) or a fraction of the chord between 0 and 1 (same as the --chord-fraction of generatePazy.py)")
    args = parser.parse_args()
    if args.twist_center not in ("axis", "auto"):
        try:
            if not 0 <= float(args.twist_center) <= 1:
                raise ValueError
        except ValueError:
            parser.error("--twist-center must be 'axis', 'auto' or a fraction of the chord between 0 and 1")

    op2_file = args.case_dir / f"{args.job}.op2"
    if not op2_file.exists():
        sys.exit(f"{op2_file} not found: run Nastran first (runNastran.py)")
    ids, disp = read_displacements(op2_file)

    meta = read_meta(args.case_dir, args.job)
    T = section_response(ids, disp, args.case_dir / meta["model"], args.stations, meta["x_axis"], meta["y_root"], meta["y_end"], args.spar_pid, args.tol, meta_sweep(meta, args.job), args.twist_center)

    out_file = args.case_dir / f"{args.job}.station_response_section{'_' + args.tag if args.tag else ''}.csv"
    np.savetxt(out_file, T.reshape(-1, 6), delimiter=",", fmt="%.17e")
    np.set_printoptions(precision=8, linewidth=150)
    print("Response at the last station:")
    print(T[-1])
    print("Saved", out_file)


if __name__ == "__main__":
    main()
