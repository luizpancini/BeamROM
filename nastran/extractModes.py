"""Extracts the eigenvalues of the modal analysis (SOL 103) from the Nastran .op2 file.

Usage:  python extractModes.py "nastran/Pazy_Technion/S0" --job S0_modal [--out S0_modal_frequencies_flextip.csv]

Writes --out (in the case directory; default <job>.eigenvalues.csv) with columns = [mode number, eigenvalue (rad/s)^2, frequency Hz] and one row per mode.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from pyNastran.op2.op2 import OP2


def modal_eigenvalues(op2_file):
    """Returns the (nmodes x 3) array with [mode number, eigenvalue (rad/s)^2, frequency Hz] of the first subcase with eigenvectors."""
    op2 = OP2(debug=False)
    op2.read_op2(str(op2_file))
    if not op2.eigenvectors:
        raise ValueError(f"{op2_file} has no eigenvectors: is it from a modal analysis (SOL 103) with DISPLACEMENT output?")
    vectors = op2.eigenvectors[min(op2.eigenvectors)]
    eigenvalues = np.asarray(vectors.eigns, dtype=float)
    modes = np.asarray(vectors.modes, dtype=float)
    return np.column_stack([modes, eigenvalues, np.sqrt(np.abs(eigenvalues)) / (2 * np.pi)])


def main():
    parser = argparse.ArgumentParser(description="Extract the eigenvalues from the .op2 of a modal analysis")
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--job", required=True, help="name of the modal job, without extension (e.g. solid_beam_modal)")
    parser.add_argument("--out", help="name of the .csv, in the case directory (default <job>.eigenvalues.csv)")
    args = parser.parse_args()

    op2_file = args.case_dir / f"{args.job}.op2"
    if not op2_file.exists():
        sys.exit(f"{op2_file} not found: run Nastran first (runNastran.py)")

    table = modal_eigenvalues(op2_file)
    out_file = args.case_dir / (args.out or f"{args.job}.eigenvalues.csv")
    np.savetxt(out_file, table, delimiter=",", fmt=["%d", "%.17e", "%.17e"])
    np.set_printoptions(precision=6, linewidth=150, suppress=True)
    print(table)
    print("Saved", out_file)


if __name__ == "__main__":
    main()
