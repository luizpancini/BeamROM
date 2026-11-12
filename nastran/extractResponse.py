"""Extracts the static response of the solid beam from the Nastran .op2 file.

Usage:  python extractResponse.py nastran/simple_beams/solid_beam/mesh1 [--job S0_static]

Writes <example>_static.response.csv (e.g. solid_beam_static.response.csv): the (6 x 6) matrix with rows = [Tx, Ty, Tz, Rx, Ry, Rz] at the load node (m and rad) and columns = load cases (Fx, Fy, Fz, Mx, My, Mz).
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from pyNastran.op2.op2 import OP2

from common import read_meta


def static_response(op2_file, load_node):
    """Returns the (6 x 6) matrix with rows = load cases (Fx, Fy, Fz, Mx, My, Mz) and columns = [Tx, Ty, Tz, Rx, Ry, Rz] at the load node (m and rad)."""
    op2 = OP2(debug=False)
    op2.read_op2(str(op2_file))
    resp = np.zeros((6, 6))
    for sid in range(1, 6 + 1):
        d = op2.displacements[sid]
        ids = d.node_gridtype[:, 0]
        idx = np.where(ids == load_node)[0][0]
        resp[sid - 1] = d.data[0, idx, :]
    return resp


def main():
    parser = argparse.ArgumentParser(description="Extract the static response from the .op2 of a case directory")
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--job", help="name of the static job, without extension (default: <example directory>_static)")
    args = parser.parse_args()

    # Jobs are named after the example directory, e.g. nastran/simple_beams/solid_beam/mesh1 -> solid_beam_static
    job = args.job or f"{args.case_dir.resolve().parent.name}_static"
    op2_file = args.case_dir / f"{job}.op2"
    if not op2_file.exists():
        sys.exit(f"{op2_file} not found: run Nastran first (runNastran.py)")

    T = static_response(op2_file, read_meta(args.case_dir, job)["load_node"]).T
    out_file = args.case_dir / f"{job}.response.csv"
    np.savetxt(out_file, T, delimiter=",", fmt="%.17e")
    np.set_printoptions(precision=8, linewidth=150)
    print(T)
    print("Saved", out_file)


if __name__ == "__main__":
    main()
