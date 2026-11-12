"""Ties the wingtip rod of the Technion Pazy wing models to the end section of the wing along the chord, as the S10 model does.

Usage:  python tieWingtipRod.py file.bdf [file2.bdf ...] [--rod-pid 13]

The S0 model ties the rod to the wing at its midpoint: one RBE2 from the grid at the center of the rod to all the grids of the end section of the wing. The rod passes through the whole chord (the S10 model glues it to the wing there), so each grid of the rod in the chord is tied by an RBE2 of its own, as in the S10 model: the independent grid is the grid of the rod and the dependent grids are those of the end section that are closest to it in the chordwise direction (x), a grid of the section to each RBE2 (all the DOF, 123456). The rod grids with no grid of the section closest to them (those beyond the chord) are not tied. The RBE2 of the rod center is replaced by those, with new ids that follow the highest id of the elements and rigid elements of the model, in the same place of the file, and the rest of the file is not changed.

The models are edited in place. A model with an RBE2 for more than one grid of the rod (e.g. S10) is left as it is, so that the script can be repeated.
"""
import argparse
from pathlib import Path

import numpy as np

from prepareSweptPazy import card_name, end_of_card, is_continuation, small_field
from tieStringerRoll import read_lines

ROD_PID = 13


def tie_wingtip_rod(lines, rod_pid=ROD_PID):
    """Lines of the bulk data with the RBE2 of the center of the rod replaced by one RBE2 for each grid of the rod in the chord (unchanged if the rod is already tied by more than one RBE2). Returns the lines and the number of RBE2 added."""
    model = read_lines(lines)
    rod = sorted({n for e in model.elements.values() if e.type == "CBEAM" and e.pid == rod_pid for n in e.node_ids})
    if not rod:
        raise SystemExit(f"No CBEAM elements with property {rod_pid} in the model")
    ties = [r for r in model.rigid_elements.values() if r.type == "RBE2" and set(r.independent_nodes) & set(rod)]
    if len(ties) != 1:
        return list(lines), 0
    tie = ties[0]
    if len(tie.independent_nodes) != 1 or str(tie.cm) != "123456":
        raise SystemExit(f"RBE2 {tie.eid} of the rod must have one independent grid and tie all the DOF")

    # Each grid of the end section goes to the grid of the rod closest to it in x
    section = sorted(tie.Gmi_node_ids)
    rod_x = np.array([model.nodes[n].xyz[0] for n in rod])
    tied = {}
    for n in section:
        tied.setdefault(rod[int(np.argmin(abs(rod_x - model.nodes[n].xyz[0])))], []).append(n)
    first = min(tied, key=lambda n: model.nodes[n].xyz[0])
    eid = max(max(model.rigid_elements, default=0), max(model.elements)) + 1
    cards = [f"$ RBE2 {eid}-{eid + len(tied) - 1}: the rod ties to the end section along the chord (tieWingtipRod.py), one for each grid of the rod in the chord"]
    for k, master in enumerate(sorted(tied, key=lambda n: model.nodes[n].xyz[0])):
        slaves = tied[master]
        cards.append(f"RBE2    {eid + k:8d}{master:8d}{'123456':>8}" + "".join(f"{n:8d}" for n in slaves[:5]))
        for i in range(5, len(slaves), 8):
            cards.append(" " * 8 + "".join(f"{n:8d}" for n in slaves[i:i + 8]))

    # Replace the card of the old tie, in the same place
    start = next(k for k, line in enumerate(lines) if not is_continuation(line) and card_name(line) == "RBE2" and "," not in line and small_field(line, 1).strip() == str(tie.eid))
    return list(lines[:start]) + cards + list(lines[end_of_card(lines, start):]), len(tied)


def main():
    parser = argparse.ArgumentParser(description="Tie the wingtip rod of the Pazy wing to the end section along the chord")
    parser.add_argument("bdf", nargs="+", type=Path)
    parser.add_argument("--rod-pid", type=int, default=ROD_PID, help="property of the CBEAM elements of the rod")
    args = parser.parse_args()
    for path in args.bdf:
        text = path.read_text()
        lines = text.split("\n")
        new, count = tie_wingtip_rod(lines, args.rod_pid)
        if not count:
            print(f"{path}: the rod is already tied by more than one RBE2")
            continue
        path.write_text("\n".join(new))
        print(f"{path}: {count} RBE2 added for the rod grids in the chord")


if __name__ == "__main__":
    main()
