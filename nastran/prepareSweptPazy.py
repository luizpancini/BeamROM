"""Prepares the bulk data of the swept Pazy wings (Technion models S10 and S20, with 10 and 20 degrees of sweep, versions LE, TE and Baseline) for the Nastran workflow of generatePazy.py.

Usage:  python prepareSweptPazy.py --sweeps 10 [20 ...] [--versions LE TE Baseline]

The sweep angles [deg] are required (there is no default): the wing of --sweeps S is in nastran/Pazy_Technion/S<S>. For each version, reads nastran/Pazy_Technion/S<S>/<version>/S<S>_<version>.nas (Femap deck, with the flutter analysis in S10, SOL SEFLUTTR, and the modal analysis in S20, SOL SEMODES) and writes nastran/Pazy_Technion/S<S>/<version>/S<S>_<version>.bdf with the bulk data only (the same layout of nastran/Pazy_Technion/S0/S0.bdf: no SPC, loads or analysis, which generatePazy.py adds). The .nas files are not modified:
  - removed: the aerodynamic model, if any (AERO, PAERO1, CAERO1, SPLINE1 and its SET1, MKAERO1/2), the flutter and eigenvalue data (FLUTTER, FLFACT, EIGRL, EIGC, PARAM LMODES, ...), the root constraint (SPC1: the root grids are clamped by generatePazy.py, which clamps the independent grid of the RBE2 of the root) and the unused Femap coordinate systems (CORD2C, CORD2S, whose ids would clash with the beam frame of generatePazy.py)
  - kept: the structural model (grids, CQUAD4, CTRIA3, CBEAM, CONM2, RBE2, properties and materials) and PARAM AUTOSPC and GRDPNT of the original deck
  - added: the RBE2 that tie the roll of the leading and trailing edge stringers to the spar at the ribs (tieStringerRoll.py, which explains them)
  - unchanged: the materials (the tip rod, PBEAM 13, and the breaking beam, PBEAM 24, are flexible). The rigid nylon of the tip (Young's modulus 1e15 Pa, with the Poisson's ratio and density of the Nylon 12) is set in the decks of the analyses by generatePazy.py --tip rigid --rigid-pids 13 24, which names the jobs with the suffix _rigidtip, and with --tip flexible they are _flextip, with these materials. A Baseline deck prepared in place by an earlier version of this script with the rigid nylon gets its original material back

The Baseline is the wing with the tip rod, without the tip weights, the fixed adapter and the calibration mass. Its deck is, for each wing:
  - its own bulk deck nastran/Pazy_Technion/S<S>/Baseline/S<S>_baseline.bdf (with the tip rod, as S10 has), if it exists and has no header (the one that this script writes), which is prepared in place, with the same removals and ENDDATA;
  - otherwise, made from the bulk data of S<S>_LE.nas (as S20, which has no Baseline deck), as the S10 Baseline differs from the S10 LE: without the tip weights, the adapter and the calibration mass (the CONM2 FPS_weights, Fixed_adapter and calib_CONM), the breaking beam (CBEAM and PBEAM 24) and its grid (the S10 Baseline keeps its calibration mass, which the S20 versions, with different ones, do not share).

The versions LE and TE differ in the tip weights and the breaking beam, and each has its own directory, since the decks and the meta.json of generatePazy.py are per directory. A file is only rewritten if its content changes, not counting its header comment (so that the Nastran results that come from the decks generatePazy.py makes from it are kept: they are removed when a deck changes).
"""
import argparse
from pathlib import Path

from common import NASTRAN_DIR
from tieStringerRoll import tie_stringer_roll

PAZY_DIR = NASTRAN_DIR / "Pazy_Technion"
HEADER_START = "$ Bulk data of "  # first line of the decks that this script writes (the own deck of a Baseline has no header)
VERSIONS = ("LE", "TE", "Baseline")
NYLON_MID = 1
RIGID_MID = 100005  # rigid nylon of --tip rigid in generatePazy.py (which an earlier version of this script set here)
BREAKING_BEAM_PID = 24
TIP_MASS_NAMES = ("FPS_weights", "Fixed_adapter", "calib_CONM")  # Femap names of the CONM2 that the Baseline made from the LE deck does not have
REMOVED_CARDS = {"AERO", "PAERO1", "CAERO1", "SPLINE1", "SET1", "MKAERO1", "MKAERO2", "FLUTTER", "FLFACT", "EIGRL", "EIGC", "SPC1", "CORD2C", "CORD2S"}
KEPT_PARAMS = {"AUTOSPC", "GRDPNT"}


def card_name(line):
    head = line.split(",")[0] if "," in line else line[:8]
    return head.strip().rstrip("*").upper()


def is_continuation(line):
    return line[:1] in ("+", "*") or not line[:8].strip()


def bulk_data(text):
    """Bulk data lines between BEGIN BULK and ENDDATA (the whole deck if it has neither, as the Baseline deck)."""
    lines = text.splitlines()
    begin = next((k + 1 for k, line in enumerate(lines) if line.strip().upper().startswith("BEGIN BULK")), 0)
    end = next((k for k, line in enumerate(lines) if line.strip().upper().startswith("ENDDATA")), len(lines))
    return lines[begin:end]


def remove_cards(lines):
    """Lines without the removed cards (and their continuation lines and the comments that precede them)."""
    out, comments, keep = [], [], True
    for line in lines:
        if line.startswith("$"):
            comments.append(line)
        elif not line.strip():
            continue
        elif is_continuation(line):
            if keep:
                out.append(line)
        else:
            name = card_name(line)
            keep = name not in REMOVED_CARDS and (name != "PARAM" or line.split(",")[1].strip().upper() in KEPT_PARAMS)
            if keep:
                out += comments + [line]
            comments = []
    return out


def small_field(line, k):
    """Field k (1 = first data field) of a line of a small-field (8 columns) card."""
    return line[8 * k:8 * (k + 1)]


def find_card(lines, name, first_field):
    """Index of the first line of the card `name` whose first data field is the integer `first_field`, in small-field format."""
    for k, line in enumerate(lines):
        if not is_continuation(line) and card_name(line) == name and "," not in line and "*" not in line[:8] and small_field(line, 1).strip() == str(first_field):
            return k
    raise SystemExit(f"{name} {first_field} (small-field format) not found")


def end_of_card(lines, k):
    k += 1
    while k < len(lines) and is_continuation(lines[k]):
        k += 1
    return k


def original_materials(lines, nylon_mid):
    """Lines without the rigid nylon of a Baseline deck that was prepared in place with it (the material RIGID_MID, an earlier version of this script): its PBEAM go back to the material `nylon_mid` (Nylon 12, which the property 13 of the model had) and the rigid material is removed, so that the deck has the original materials (the rigid ones are set by generatePazy.py with --tip rigid)."""
    mats = [k for k, line in enumerate(lines) if not is_continuation(line) and card_name(line) == "MAT1" and small_field(line, 1).strip() == str(RIGID_MID)]
    if not mats:
        return lines
    lines = list(lines)
    for k, line in enumerate(lines):
        if not is_continuation(line) and card_name(line) == "PBEAM" and small_field(line, 2).strip() == str(RIGID_MID):
            lines[k] = line[:16] + f"{nylon_mid:>8d}" + line[24:]
    k = mats[0]
    del lines[k:end_of_card(lines, k)]
    if k > 0 and lines[k - 1].startswith(f"$ Material {RIGID_MID}"):
        del lines[k - 1]
    return lines


def has_header(bdf):
    return bdf.exists() and bdf.read_text().startswith(HEADER_START)


def version_files(sweep, version, source_dir):
    """Source deck and prepared bulk data of a version of the wing with `sweep` degrees of sweep, in its directory: the Femap deck (S<sweep>_<version>.nas) of LE and TE, and for the Baseline its own bulk deck, prepared in place, or else the LE deck, from which it is made."""
    if version == "Baseline":
        bdf = source_dir / version / f"S{sweep}_baseline.bdf"
        return (bdf if bdf.exists() and not has_header(bdf) else source_dir / "LE" / f"S{sweep}_LE.nas"), bdf
    return source_dir / version / f"S{sweep}_{version}.nas", source_dir / version / f"S{sweep}_{version}.bdf"


def card_id(line):
    return small_field(line, 1).strip()


def baseline_from_tip_version(lines):
    """Lines of the bulk data of the Baseline made from those of a version with the tip masses: without the CONM2 of the tip weights, of the adapter and of the calibration mass (named in the Femap comment before them), the CBEAM of the breaking beam, its property and the grids that only it uses (the end that is not in the rod)."""
    cbeams = [line for line in lines if not is_continuation(line) and card_name(line) == "CBEAM"]
    breaking = [line for line in cbeams if small_field(line, 2).strip() == str(BREAKING_BEAM_PID)]
    rod = {small_field(line, k).strip() for line in cbeams if line not in breaking for k in (3, 4)}
    free_ends = {small_field(line, k).strip() for line in breaking for k in (3, 4)} - rod
    out, comments, keep = [], [], True
    for line in lines:
        if line.startswith("$"):
            comments.append(line)
            continue
        if is_continuation(line):
            if keep:
                out.append(line)
            continue
        name = card_name(line)
        label = comments[-1] if comments else ""
        drop = (name == "CONM2" and any(w in label for w in TIP_MASS_NAMES)) \
            or (name == "CBEAM" and small_field(line, 2).strip() == str(BREAKING_BEAM_PID)) \
            or (name == "PBEAM" and card_id(line) == str(BREAKING_BEAM_PID)) \
            or (name == "GRID" and card_id(line) in free_ends)
        keep = not drop
        if keep:
            out += comments + [line]
        comments = []
    return out + comments


def prepare(sweep, version, source_dir, nylon_mid):
    nas, bdf = version_files(sweep, version, source_dir)
    if not nas.exists():
        raise SystemExit(f"{nas} not found")
    lines = remove_cards(bulk_data(nas.read_text()))
    made = nas != bdf and version == "Baseline"  # Baseline made from the LE deck (it has no deck of its own)
    if made:
        lines = baseline_from_tip_version(lines)
    lines = original_materials(lines, nylon_mid)
    # Roll of the stringers tied to the spar at the ribs (already in the Baseline deck, if prepared in place)
    lines, _ = tie_stringer_roll(lines)
    if nas == bdf:
        text = "\n".join(lines + ["ENDDATA", ""])
    else:
        note = "; tip weights, adapter, calibration mass and breaking beam removed" if made else ""
        header = [f"$ Bulk data of {nas.name} prepared by prepareSweptPazy.py (aerodynamic and flutter data, eigenvalue methods and constraints removed{note})"]
        text = "\n".join(header + lines + ["ENDDATA", ""])
    bdf.parent.mkdir(parents=True, exist_ok=True)
    # Not rewritten if only the header comment changes (e.g. the name of this script, which was preparePazyS10.py)
    old = bdf.read_text().split("\n")[1:] if bdf.exists() and nas != bdf else None
    if not bdf.exists() or (old != text.split("\n")[1:] if old is not None else bdf.read_text() != text):
        bdf.write_text(text)
        print("Written", bdf)
    return bdf


def main():
    parser = argparse.ArgumentParser(description="Prepare the bulk data of the swept Pazy wings (S10, S20) for the Nastran workflow")
    parser.add_argument("--sweeps", nargs="+", type=int, required=True, help="sweep angles [deg] of the wings (integers: the wing is in nastran/Pazy_Technion/S<sweep>)")
    parser.add_argument("--versions", nargs="+", choices=VERSIONS, default=list(VERSIONS))
    parser.add_argument("--source-dir", type=Path, default=None, help="directory with a subdirectory per version, each with its S<sweep>_<version>.nas file (default: nastran/Pazy_Technion/S<sweep>)")
    parser.add_argument("--nylon-mid", type=int, default=NYLON_MID, help="material (Nylon 12) of the property 13 (tip rod), restored in a Baseline deck prepared in place with the rigid nylon")
    args = parser.parse_args()
    for sweep in args.sweeps:
        source_dir = (args.source_dir or PAZY_DIR / f"S{sweep}").resolve()
        for version in args.versions:
            prepare(sweep, version, source_dir, args.nylon_mid)


if __name__ == "__main__":
    main()
