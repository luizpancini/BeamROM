"""Extracts the generalized displacements of the clamped windows of the discretized beam of the Pazy wing from the Nastran .op2 of the segments deck (generatePazy.py --analyses segments).

Usage:  python extractSegments.py "nastran/Pazy_Technion/S0" --job S0_segments

In the deck, subcases 6(w-1)+1 ... 6w are the six unit loads [F1, F2, F3, M1, M2, M3] at the grid of the station w+m, with the section of the station w clamped, where m is the number of intervals per window (`intervals` in <job>.json). The displacements of that grid are those of the end of the window w, a cantilever clamped at its upstream node and loaded at its downstream node. The stations come from <job>.json, written by generatePazy.py.

Writes <job>.segment_response.csv: the windows stacked, (6 (n-m) x 6) matrix where rows 6(w-1)+1:6w are [Tx, Ty, Tz, Rx, Ry, Rz] of the downstream node of window w (m and rad, in the beam frame: x1 = y of the Nastran model along the span, x2 = -x, x3 = z) and columns are the load cases (Fx, Fy, Fz, Mx, My, Mz).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from pyNastran.op2.op2 import OP2


def read_displacements(op2_file, n_subcases):
    """Displacements of the subcases: (grid ids, array (n_subcases, n grids, 6))."""
    op2 = OP2(debug=False)
    op2.read_op2(str(op2_file))
    missing = [sid for sid in range(1, n_subcases + 1) if sid not in op2.displacements]
    if missing:
        sys.exit(f"Subcases {missing} not found in {op2_file}")
    ids = op2.displacements[1].node_gridtype[:, 0]
    return ids, np.array([op2.displacements[sid].data[0] for sid in range(1, n_subcases + 1)], dtype=float)


def segment_response(ids, disp, stations, intervals=1):
    """(n-m, 6, 6) generalized displacements of the downstream node of each window of `intervals` = m intervals, already in the beam frame (displacement coordinate system of the grids of the stations)."""
    rows = {int(g): i for i, g in enumerate(ids)}
    n_windows = len(stations) - intervals
    T = np.zeros((n_windows, 6, 6))
    for k in range(n_windows):
        T[k] = disp[6 * k:6 * (k + 1), rows[stations[k + intervals]["grid"]], :].T
    return T


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--job", required=True, help="name of the segments job, without extension (e.g. S0_segments)")
    args = parser.parse_args()

    op2_file = args.case_dir / f"{args.job}.op2"
    if not op2_file.exists():
        sys.exit(f"{op2_file} not found: run Nastran first (runNastran.py)")
    stations_file = args.case_dir / f"{args.job}.json"
    if not stations_file.exists():
        sys.exit(f"{stations_file} not found: generate the segments deck first (generatePazy.py --analyses segments)")
    meta = json.loads(stations_file.read_text())
    stations, intervals = meta["stations"], meta.get("intervals", 1)

    ids, disp = read_displacements(op2_file, 6 * (len(stations) - intervals))
    T = segment_response(ids, disp, stations, intervals)

    out_file = args.case_dir / f"{args.job}.segment_response.csv"
    np.savetxt(out_file, T.reshape(-1, 6), delimiter=",", fmt="%.17e")
    np.set_printoptions(precision=8, linewidth=150)
    print("Response of the last window:")
    print(T[-1])
    print("Saved", out_file)


if __name__ == "__main__":
    main()
