__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

# M04. Lever arm — no dedicated QTCAD tutorial script exists for this topic.
# This script applies qtcad.device.leverarm / leverarm_matrix directly to the
# M03 double-dot FD-SOI device (get_double_dot_fdsoi), following the usage
# pattern shown in examples/practical_application/GaAs_gated/4-lever_arm.py.

import pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from qtcad.device import constants as ct
from qtcad.device.mesh3d import Mesh
from qtcad.device.poisson import SolverParams as PoissonSolverParams
from qtcad.device.leverarm import Solver as LeverArmSolver
from qtcad.device.leverarm import SolverParams as LeverArmSolverParams
from qtcad.device.leverarm_matrix import Solver as LeverArmMatrixSolver
from qtcad.device.leverarm_matrix import SolverParams as LeverArmMatrixSolverParams
from helper.double_dot_fdsoi import get_double_dot_fdsoi

script_dir = pathlib.Path(__file__).parent.resolve()
out_dir = script_dir / "output"
out_dir.mkdir(exist_ok=True)

# NOTE on mesh choice: this script needs ~11 independent Poisson+Schrodinger
# solves (5 for the plunger lever-arm sweep + 6 for the 5-gate finite-difference
# cross-coupling matrix), none of which can reuse a warm-started potential the
# way M03/M16 do. Using the same *refined* mesh as M03 (45 MB, adaptively
# refined for a single production point) would make this prohibitively slow.
# The lever arm is a *slope* (dE/dV), which is far less sensitive to mesh
# refinement than the absolute confinement energy M03 reports, so the
# un-refined base mesh (dqdfdsoi.msh, same geometry/physical regions, no
# adaptive refinement) is used here instead. This is a deliberate accuracy/
# speed trade-off for this module only, not a change to M03's own mesh.
path_mesh = script_dir / "meshes" / "dqdfdsoi.msh"
scaling = 1e-9
mesh = Mesh(scaling, str(path_mesh))

# Reference bias point: same "low barrier" tuning M03/M16 use (no detuning).
back_gate_bias = -0.5
barrier_gate_1_bias = 0.5
plunger_gate_1_bias = 0.59
barrier_gate_2_bias = 0.57
plunger_gate_2_bias = 0.59
barrier_gate_3_bias = 0.5

dot_region_list = ["oxide_dot", "gate_oxide_dot", "buried_oxide_dot", "channel_dot"]

poisson_params = PoissonSolverParams({"tol": 1e-3})

# -----------------------------------------------------------------------------
# 1. Single-gate lever arm of plunger_gate_1 (polynomial/linear fit of E0 vs V)
# -----------------------------------------------------------------------------
dvc1 = get_double_dot_fdsoi(
    mesh, back_gate_bias, barrier_gate_1_bias, plunger_gate_1_bias,
    barrier_gate_2_bias, plunger_gate_2_bias, barrier_gate_3_bias,
)

Vp_vals = np.linspace(0.56, 0.62, 5)  # sweep around the reference bias (0.59 V)
leverarm_params = LeverArmSolverParams({"pot_solver_params": poisson_params})

print("Computing lever arm of plunger_gate_1_bnd ...")
leverarm_slv = LeverArmSolver(
    dvc1,
    "plunger_gate_1_bnd",
    Vp_vals,
    dot_region=dot_region_list,
    solver_params=leverarm_params,
)
poly_coeffs = leverarm_slv.solve(state=0, degree=1)
lever_arm_p1 = abs(poly_coeffs[0]) / ct.e  # dimensionless (eV/V)
print(f"Lever arm of plunger_gate_1 on E0: {lever_arm_p1:.4f}")

leverarm_slv.save(str(out_dir / "leverarm_plunger_gate_1.txt"))

fig, ax = plt.subplots()
ax.set_xlabel("plunger_gate_1 voltage (V)")
ax.set_ylabel("Energy (eV)")
for i, data in enumerate(leverarm_slv.energies.T / ct.e):
    ax.plot(Vp_vals, data, "o-", label=f"state {i}")
ax.legend()
fig.tight_layout()
fig.savefig(str(out_dir / "leverarm_plunger_gate_1.png"), dpi=150)

# -----------------------------------------------------------------------------
# 2. Cross-coupling lever-arm matrix across all 5 gates (finite difference)
# -----------------------------------------------------------------------------
dvc2 = get_double_dot_fdsoi(
    mesh, back_gate_bias, barrier_gate_1_bias, plunger_gate_1_bias,
    barrier_gate_2_bias, plunger_gate_2_bias, barrier_gate_3_bias,
)

labels = [
    "barrier_gate_1_bnd", "plunger_gate_1_bnd", "barrier_gate_2_bnd",
    "plunger_gate_2_bnd", "barrier_gate_3_bnd",
]
ref_potentials = [
    barrier_gate_1_bias, plunger_gate_1_bias, barrier_gate_2_bias,
    plunger_gate_2_bias, barrier_gate_3_bias,
]

matrix_params = LeverArmMatrixSolverParams({"pot_solver_params": poisson_params})

print("\nComputing lever arm matrix (5 gates, finite difference) ...")
matrix_slv = LeverArmMatrixSolver(
    dvc2, labels, ref_potentials, dot_region=dot_region_list,
    solver_params=matrix_params,
)
lever_arm_matrix = matrix_slv.solve(bias_increment=1e-3)
lever_arm_matrix_eV_per_V = np.abs(lever_arm_matrix) / ct.e

print("Lever arm matrix (rows: eigenstates 0,1 ; cols: gates):")
print(labels)
print(lever_arm_matrix_eV_per_V[:2])

matrix_slv.save(str(out_dir / "leverarm_matrix.txt"))
np.savetxt(str(out_dir / "leverarm_matrix_eV_per_V.txt"), lever_arm_matrix_eV_per_V)
with open(out_dir / "leverarm_matrix_labels.txt", "w") as f:
    f.write("\n".join(labels))

print("\nDone.")
