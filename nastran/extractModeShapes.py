"""Extracts the mode shapes of the Technion Pazy wing (straight, S0, or swept, S10 and S20) from the Nastran .op2 of the modal analysis (SOL 103), to compare them with those of the beam model (modal assurance criterion).

Usage:  python extractModeShapes.py "nastran/Pazy_Technion/S0" --job S0_modal_flextip [--nmodes 15] [--name S0_modal_shapes_flextip]   (swept: "nastran/Pazy_Technion/S10/LE" --job S10_LE_modal_flextip --name S10_LE_modal_shapes_flextip)

Writes <name>_wing.csv and <name>_rod.csv (--name in the case directory; default <job>.shapes), with the translations (basic frame of the model: x = chord, positive aft, y = span, z = vertical) of the first --nmodes modes, each normalized so that its largest translation (of the points of both files) has magnitude 1:
  <name>_wing.csv  columns mode, y [m, from the root section], uLE_x, uLE_y, uLE_z, uTE_x, uTE_y, uTE_z, one row per row of spar grids (those at the same y, within --tol, with at least MIN_SECTION_GRIDS grids) of each mode. The translations are those of the leading edge (LE) and trailing edge (TE) points of the section (at the minimum and maximum x of the shell grids of the model, on the plane z = 0; in a swept wing, shifted along x by y tan(sweep) with the sections, with y the one of the basic frame of the model), from the rigid body motion of the section fitted to the translations of the spar grids (`section_motion` of extractStations.py). The rotation of the section about the chord axis does not move those points
  <name>_rod.csv   columns mode, x [m, from the leading edge], u_x, u_y, u_z, one row per grid of the tip rod (CBEAM of property --rod-pid) of each mode
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2

from extractStations import MIN_SECTION_GRIDS, SPAR_PID, TOL, section_motion, spar_grids, spar_rows
from generatePazy import ROD_PID, beam_axis_x


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--job", required=True, help="name of the modal job, without extension (e.g. S0_modal_flextip)")
    parser.add_argument("--nmodes", type=int, default=15, help="number of modes to extract")
    parser.add_argument("--name", help="name of the outputs, without the _wing.csv and _rod.csv, in the case directory (default <job>.shapes)")
    parser.add_argument("--spar-pid", type=int, default=SPAR_PID, help="property of the shell elements of the spar")
    parser.add_argument("--rod-pid", type=int, default=ROD_PID, help="property of the CBEAM elements of the tip rod")
    parser.add_argument("--tol", type=float, default=TOL, help="tolerance in y to take the grids of the spar at a section")
    args = parser.parse_args()

    op2_file = args.case_dir / f"{args.job}.op2"
    if not op2_file.exists():
        sys.exit(f"{op2_file} not found: run Nastran first (runNastran.py)")
    bdf_file = args.case_dir / f"{args.job}.bdf"
    if not bdf_file.exists():
        sys.exit(f"{bdf_file} not found: generate the modal deck first (generatePazy.py --analyses modal)")

    model = BDF(debug=False)
    model.read_bdf(str(bdf_file), xref=False, punch=False)
    op2 = OP2(debug=False)
    op2.read_op2(str(op2_file))
    if not op2.eigenvectors:
        sys.exit(f"{op2_file} has no eigenvectors: is it from a modal analysis (SOL 103) with DISPLACEMENT output?")
    vectors = op2.eigenvectors[min(op2.eigenvectors)]
    nmodes = min(args.nmodes, vectors.data.shape[0])
    rows = {int(g): i for i, g in enumerate(vectors.node_gridtype[:, 0])}
    disp = np.asarray(vectors.data, dtype=float)[:nmodes]  # (modes, grids, 6)
    if any(model.nodes[n].cd != 0 for n in rows if n in model.nodes):
        sys.exit("Grids with a non-basic displacement coordinate system are not supported")

    # Sections of the spar, with the leading and trailing edge points on the plane z = 0
    spar = spar_grids(model, args.spar_pid)
    sweep = json.loads(bdf_file.with_suffix(".json").read_text())["sweep"]  # written by generatePazy.py next to the modal deck
    x_le, x_te = beam_axis_x(model, 0.0, sweep), beam_axis_x(model, 1.0, sweep)
    tan_sweep = np.tan(np.radians(sweep))
    y_root = min(model.nodes[n].xyz[1] for n in spar)
    section_y, sections = [], []
    for y in spar_rows(model, spar, args.tol):
        section = [n for n in spar if abs(model.nodes[n].xyz[1] - y) < args.tol]
        if len(section) < MIN_SECTION_GRIDS:
            continue
        points = np.array([model.nodes[n].xyz for n in section])
        idx = [rows[n] for n in section]
        x_le_y = x_le + y * tan_sweep  # the sections of a swept wing are shifted along x
        p_le = np.array([x_le_y, y, 0.])
        # Translation of the LE point, and of the TE point (which is at the distance x_te - x_le along x from it)
        u = np.array([section_motion(points, disp[k, idx, :], p_le, x_le_y) for k in range(nmodes)])  # (modes, [u, θ])
        r = np.array([x_te - x_le, 0., 0.])
        u_te = u[:, :3] + np.cross(u[:, 3:], r)
        section_y.append(y - y_root)
        sections.append(np.hstack([u[:, :3], u_te]))
    sections = np.array(sections).transpose(1, 0, 2)  # (modes, sections, 6)

    # Grids of the rod, along x
    rod = sorted({n for e in model.elements.values() if e.type == "CBEAM" and e.pid == args.rod_pid for n in e.node_ids}, key=lambda n: model.nodes[n].xyz[0])
    if not rod:
        sys.exit(f"No CBEAM elements with property {args.rod_pid} in the model")
    rod_x = np.array([model.nodes[n].xyz[0] - x_le - tan_sweep * model.nodes[n].xyz[1] for n in rod])
    rod_u = disp[:, [rows[n] for n in rod], :3]  # (modes, grids, 3)

    # Normalization of each mode: the largest translation of the points of the files has magnitude 1
    scale = np.maximum(np.abs(sections).max(axis=(1, 2)), np.abs(rod_u).max(axis=(1, 2)))
    sections /= scale[:, None, None]
    rod_u /= scale[:, None, None]

    wing_rows = [[k + 1, y, *sections[k, i]] for k in range(nmodes) for i, y in enumerate(section_y)]
    rod_rows = [[k + 1, x, *rod_u[k, i]] for k in range(nmodes) for i, x in enumerate(rod_x)]
    stem = args.name or f"{args.job}.shapes"
    for kind, header, table, fmt in (
        ("wing", "mode,y,uLE_x,uLE_y,uLE_z,uTE_x,uTE_y,uTE_z", wing_rows, ["%d"] + ["%.10e"] * 7),
        ("rod", "mode,x,u_x,u_y,u_z", rod_rows, ["%d"] + ["%.10e"] * 4),
    ):
        out_file = args.case_dir / f"{stem}_{kind}.csv"
        np.savetxt(out_file, np.array(table), delimiter=",", fmt=fmt, header=header, comments="")
        print("Saved", out_file)


if __name__ == "__main__":
    main()
