"""Runs Nastran on every .bdf of a case directory (standard library only, so it also runs on a Windows machine/VM that has Nastran but not the project's Python environment).

Usage:  python runNastran.py nastran/simple_beams/solid_beam/mesh1 [--jobs solid_beam_static] [--force]

Jobs whose .op2 already exists are skipped. If the Nastran executable is not found (set NASTRAN_EXE or --exe, or put it on PATH), the pending .bdf files are listed and the script exits with code 2: run them on a machine with Nastran, copy the outputs (.op2) back to the same directory and continue the workflow.
"""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

EXIT_NO_SOLVER = 2


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--exe", default=os.environ.get("NASTRAN_EXE", "nastran"), help="Nastran executable (default: $NASTRAN_EXE or 'nastran')")
    parser.add_argument("--nastran-args", nargs="*", default=["old=no", "scr=yes", "batch=no"], help="extra arguments passed to Nastran")
    parser.add_argument("--force", action="store_true", help="rerun even if the .op2 already exists")
    parser.add_argument("--jobs", nargs="*", help="names of the .bdf files to run, without extension (default: all)")
    args = parser.parse_args()

    case_dir = args.case_dir.resolve()
    bdfs = sorted(bdf for bdf in case_dir.glob("*.bdf") if not args.jobs or bdf.stem in args.jobs)
    if not bdfs:
        sys.exit(f"No matching .bdf files in {case_dir}")

    pending = [bdf for bdf in bdfs if args.force or not bdf.with_suffix(".op2").exists()]
    if not pending:
        print("All .op2 files already exist, nothing to run.")
        return

    exe = shutil.which(args.exe)
    if exe is None:
        print(f"Nastran executable '{args.exe}' not found. Run these files on a machine with Nastran:", file=sys.stderr)
        for bdf in pending:
            print(f"  {bdf}", file=sys.stderr)
        print(f"and copy the resulting .op2 files to {case_dir}", file=sys.stderr)
        sys.exit(EXIT_NO_SOLVER)

    for bdf in pending:
        print(f"Running Nastran on {bdf.name} ...", flush=True)
        result = subprocess.run([exe, bdf.name, *args.nastran_args], cwd=case_dir)
        if result.returncode != 0 or not bdf.with_suffix(".op2").exists():
            sys.exit(f"Nastran failed on {bdf.name} (exit code {result.returncode}); check {bdf.with_suffix('.f06').name}")
    print("Nastran finished.")


if __name__ == "__main__":
    main()
