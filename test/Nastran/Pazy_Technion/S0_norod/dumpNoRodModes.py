"""Dumps, for the first modes of the Nastran modal analysis of the S0 Pazy wing without the tip rod (SOL 103, clamped at the root), the displacements of all grids and the rigid body motion of each section of the spar (as `section_motion` of BeamROM extractStations.py), both in the basic frame of the model, to compute the generalized masses of the modes with the equivalent beam mass models (modal_mass_study.jl).

Usage:  <BeamROM python> dumpNoRodModes.py <BeamROM nastran folder> <output folder> [nmodes] [reference_axis]
Writes S0_norod_modal_grids.csv (mode, grid, T1..T3, R1..R3) and S0_norod_modal_sections.csv (mode, y, u1..u3, θ1..θ3, at the reference axis, with x = reference_axis of the chord and z = 0).
"""
import sys
from pathlib import Path

import numpy as np

nastran_dir, out_dir = Path(sys.argv[1]), Path(sys.argv[2])
nmodes = int(sys.argv[3]) if len(sys.argv) > 3 else 15
reference_axis = float(sys.argv[4]) if len(sys.argv) > 4 else 0.4653
sys.path.insert(0, str(nastran_dir))

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2 import OP2

from extractStations import MIN_SECTION_GRIDS, SPAR_PID, TOL, section_motion, spar_grids, spar_rows
from generatePazy import beam_axis_x

case = nastran_dir / "Pazy_Technion" / "S0_norod"
model = BDF(debug=False)
model.read_bdf(str(case / "S0_norod_modal.bdf"), xref=False, punch=False)
op2 = OP2(debug=False)
op2.read_op2(str(case / "S0_norod_modal.op2"))
vectors = op2.eigenvectors[min(op2.eigenvectors)]
nmodes = min(nmodes, vectors.data.shape[0])
ids = [int(g) for g in vectors.node_gridtype[:, 0]]
rows = {g: i for i, g in enumerate(ids)}
disp = np.asarray(vectors.data, dtype=float)[:nmodes]  # (modes, grids, 6)

with open(out_dir / "S0_norod_modal_grids.csv", "w") as f:
    f.write("mode,grid,T1,T2,T3,R1,R2,R3\n")
    for k in range(nmodes):
        for i, g in enumerate(ids):
            f.write(f"{k+1},{g}," + ",".join(f"{v:.12e}" for v in disp[k, i]) + "\n")

spar = spar_grids(model, SPAR_PID)
x_ref = beam_axis_x(model, reference_axis, 0.0)
with open(out_dir / "S0_norod_modal_sections.csv", "w") as f:
    f.write("mode,y,u1,u2,u3,th1,th2,th3\n")
    for y in spar_rows(model, spar, TOL):
        section = [n for n in spar if abs(model.nodes[n].xyz[1] - y) < TOL]
        if len(section) < MIN_SECTION_GRIDS:
            continue
        points = np.array([model.nodes[n].xyz for n in section])
        idx = [rows[n] for n in section]
        for k in range(nmodes):
            m = section_motion(points, disp[k, idx, :], np.array([x_ref, y, 0.0]), x_ref)
            f.write(f"{k+1},{y:.10f}," + ",".join(f"{v:.12e}" for v in m) + "\n")
print("done", nmodes, "modes,", len(ids), "grids")
