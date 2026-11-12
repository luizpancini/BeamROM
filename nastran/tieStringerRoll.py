"""Ties the roll of the leading and trailing edge stringers of the Technion Pazy wing models to the spar at the ribs.

Usage:  python tieStringerRoll.py file.bdf [file2.bdf ...]

The stringers (CBEAM of properties 1, leading edge, and 2, trailing edge) are connected to the rest of the model only by the 12 micron skin and, at each rib, by one vertex of the shell of the rib (properties 4 and 10, in planes of constant y), whose rotation about the normal of the rib (the span direction y, the drilling rotation of the shell, which has almost no stiffness) is the roll of the stringer. The stringers could then roll freely relative to the sections of the wing, which gives local modes of the stringers and lowers the torsion frequencies by up to 40% relative to the equivalent beam, whose sections are rigid. The cards added here tie that rotation at each rib to the rotation of the grid of the edge of the spar (minimum x for the leading edge and maximum x for the trailing edge) at the same y: one RBE2 for each grid of the stringer shared with the shell of a rib, with only the rotation about y (component 5, in the basic frame, in which all grids are defined) dependent. The stringers can still roll between the ribs, but not at them, as if the stringers were part of the ribs. The static response does not change (the roll of the stringers follows the twist of the section), and neither do the stiffness and inertia matrices of the equivalent beam, which come from the static analysis and from the mass of the model.

The models are edited in place: the cards are added before ENDDATA (or at the end of the file, if it has none), after a comment line that is also the mark of the models that already have them, so that the script can be repeated. The ids of the RBE2 follow the highest id of the elements and rigid elements of the model. `tie_stringer_roll(lines)` does the same on the lines of the bulk data, for prepareSweptPazy.py.
"""
import sys
import tempfile
from pathlib import Path

import numpy as np
from pyNastran.bdf.bdf import BDF

SPAR_PID = 100003
STRINGER_PIDS = (1, 2)  # leading edge, trailing edge
RIB_PIDS = {1: 4, 2: 10}  # shells of the ribs that share the grids of each stringer at the ribs
ROLL_DOF = "5"  # rotation about y (the span, the normal of the ribs) in the basic frame
TOL = 1e-5
MARK = "$ Roll of the leading and trailing edge stringers tied to the spar at the ribs (tieStringerRoll.py): RBE2 of the rotation about y"


def read_lines(lines):
    """pyNastran model of the bulk data given by its lines."""
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "model.bdf"
        path.write_text("\n".join(lines) + "\n")
        model = BDF(debug=False)
        model.read_bdf(str(path), xref=False, punch=True)
    return model


def roll_ties(model, spar_pid=SPAR_PID):
    """Pairs (grid of the edge of the spar, grid of the stringer at a rib) to tie, in the order of the stringers and of y."""
    spar = sorted({n for e in model.elements.values() if e.type in ("CQUAD4", "CTRIA3") and e.pid == spar_pid for n in e.node_ids})
    if not spar:
        raise SystemExit(f"No shell elements with property {spar_pid} in the model")
    ys = [model.nodes[n].xyz[1] for n in model.nodes]
    y_root, y_end = min(ys), max(model.nodes[n].xyz[1] for n in spar)
    dependent = {n for r in model.rigid_elements.values() for n in r.Gmi_node_ids}
    pairs = []
    for pid in STRINGER_PIDS:
        beams = [e for e in model.elements.values() if e.type == "CBEAM" and e.pid == pid]
        if not beams:
            raise SystemExit(f"No CBEAM elements with property {pid} (stringer) in the model")
        stringer = {n for e in beams for n in e.node_ids}
        rib = {n for e in model.elements.values() if e.type in ("CQUAD4", "CTRIA3") and e.pid == RIB_PIDS[pid] for n in e.node_ids}
        for n in sorted(stringer & rib, key=lambda n: (model.nodes[n].xyz[1], n)):
            x, y, _ = model.nodes[n].xyz
            if abs(y - y_root) < TOL or abs(y - y_end) < TOL:
                raise SystemExit(f"Grid {n} of the stringer {pid} is at a clamped or tied end section (y = {y})")
            edge = [m for m in spar if abs(model.nodes[m].xyz[1] - y) < TOL]
            if not edge:
                raise SystemExit(f"No spar grid at y = {y} (stringer {pid}, grid {n})")
            master = (min if pid == 1 else max)(edge, key=lambda m: model.nodes[m].xyz[0])
            if n in dependent or master in dependent:
                raise SystemExit(f"Grid {n} (stringer {pid}) or {master} (spar) is dependent in a rigid element")
            pairs.append((pid, master, n))
    return pairs


def tie_stringer_roll(lines):
    """Lines of the bulk data with the RBE2 of the roll of the stringers added at the end (unchanged if it has them). Returns the lines and the number of RBE2 added."""
    if any(line.strip() == MARK for line in lines):
        return lines, 0
    model = read_lines(lines)
    pairs = roll_ties(model)
    eid = max(max(model.rigid_elements, default=0), max(model.elements)) + 1
    cards = [MARK] + [f"RBE2    {eid + k:8d}{master:8d}{ROLL_DOF:>8}{n:8d}" for k, (_, master, n) in enumerate(pairs)]
    return list(lines) + cards, len(pairs)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for name in sys.argv[1:]:
        path = Path(name)
        lines = path.read_text().split("\n")
        if lines and lines[-1] == "":
            lines.pop()
        end = next((k for k, line in enumerate(lines) if line.strip().upper().startswith("ENDDATA")), len(lines))
        new, count = tie_stringer_roll(lines[:end])
        if not count:
            print(f"{path}: already has the ties")
            continue
        path.write_text("\n".join(new + lines[end:]) + "\n")
        print(f"{path}: {count} RBE2 added")


if __name__ == "__main__":
    main()
