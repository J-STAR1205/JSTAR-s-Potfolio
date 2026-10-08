__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

from qtcad.device import materials as mt
from qtcad.device import io
from qtcad.device.mesh3d import Mesh
from qtcad.device import Device
from qtcad.device.leverarm_matrix import Solver as LeverArmSolver
from qtcad.device.leverarm_matrix import SolverParams as LeverArmSolverParams
from qtcad.device.poisson import SolverParams as PoissonSolverParams
from qtcad.device.schrodinger import SolverParams as SchrodingerSolverParams
import numpy as np
from pathlib import Path

script_dir = Path(__file__).parent.resolve()
mesh_dir = script_dir / "meshes"
mesh_file_dir = mesh_dir / "refined_gedqd.msh"
xao_dir = mesh_dir / "ge_dqd.xao"
out_dir = script_dir / "output"
out_dir.mkdir(exist_ok=True)

scaling = 1e-9  # nanometers
mesh = Mesh(scaling, mesh_file_dir)
# mesh.show()  # disabled: crashes headless (see note in 2-poisson_schrod.py)

dvc = Device(mesh, conf_carriers="h", hole_kp_model="luttinger_kohn_foreman")
dvc.set_temperature(0.1)

# Assign the material to regions
dvc.new_region("substrate", mt.SiGe_DFT)
dvc.new_region("SiGe_barrier", mt.SiGe_DFT)
dvc.new_region("Ge_well", mt.Ge)
dvc.new_region("SiGe_cap", mt.SiGe_DFT)
dvc.new_region("oxide", mt.Al2O3)

# Dot-region volumes
dvc.new_region("Ge_well.dot_region", mt.Ge)
dvc.new_region("SiGe_barrier.dot_region", mt.SiGe_DFT)
dvc.new_region("SiGe_cap.dot_region", mt.SiGe_DFT)

# boundary conditions
Ew = mt.Ge.chi + 1.1 * mt.Ge.Eg
dvc.new_gate_bnd("P1", -0.6, Ew)
dvc.new_gate_bnd("P2", -0.6, Ew)
dvc.new_gate_bnd("BC", 0.9, Ew)
dvc.new_gate_bnd("BL", 0.5, Ew)
dvc.new_gate_bnd("BR", 0.5, Ew)

dot_region_list = [
    "Ge_well.dot_region",
    "SiGe_barrier.dot_region",
    "SiGe_cap.dot_region",
]

dvc.set_dot_region(dot_region_list)

phi = io.load(out_dir / "electrostatic_potential.hdf5")
dvc.set_potential(phi)

dvc.set_V_from_phi()

bias_vector = np.array([-1.5, -1.5])
gate_labels = ["P1", "P2"]

poisson_params = PoissonSolverParams()
poisson_params.tol = 1e-3
poisson_params.maxiter = 50

schrodinger_params = SchrodingerSolverParams()
schrodinger_params.num_states = 20
schrodinger_params.tol = 1e-5  # eV

lam_params = LeverArmSolverParams()
lam_params.pot_solver_params = poisson_params
lam_params.schrod_solver_params = schrodinger_params

bias_inc = 1e-3

slv = LeverArmSolver(
    dvc=dvc,
    labels=gate_labels,
    potentials=bias_vector,
    dot_region=dot_region_list,
    solver_params=lam_params,
)

lever_arm_matrix = slv.solve(bias_inc)

print("Lever-arm matrix")
print(lever_arm_matrix)

np.save(out_dir / "lever_arm_matrix.npy", lever_arm_matrix)
