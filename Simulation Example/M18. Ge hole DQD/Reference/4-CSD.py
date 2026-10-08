__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

from qtcad.device import materials as mt
from qtcad.device import io
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device import Device, SubDevice
from qtcad.device.many_body import Solver as ManyBodySolver
from qtcad.device.many_body import SolverParams as ManyBodySolverParams
from qtcad.device.schrodinger import Solver as SchrodingerSolver
from qtcad.device.schrodinger import SolverParams as SchrodingerSolverParams
from qtcad.transport.junction import Junction
from qtcad.transport.mastereq import add_spectrum
import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt

script_dir = Path(__file__).parent.resolve()
mesh_dir = script_dir / "meshes"
mesh_file_dir = mesh_dir / "refined_gedqd.msh"
xao_dir = mesh_dir / "ge_dqd.xao"
out_dir = script_dir / "output"
out_dir.mkdir(exist_ok=True)

scaling = 1e-9  # nanometers
mesh = Mesh(scaling, mesh_file_dir)

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

submesh = SubMesh(dvc.mesh, dot_region_list)
subdvc = SubDevice(dvc, submesh)

schrodinger_params = SchrodingerSolverParams()
schrodinger_params.num_states = 20
schrodinger_params.tol = 1e-5

schrodinger_solver = SchrodingerSolver(subdvc, solver_params=schrodinger_params)
schrodinger_solver.solve()

subdvc.print_energies()
energies = subdvc.energies

many_body_solver_params = ManyBodySolverParams()
many_body_solver_params.n_degen = 1
many_body_solver_params.num_states = 8

slv = ManyBodySolver(subdvc, solver_params=many_body_solver_params)

coulomb_overlap = slv.get_coulomb_matrix(overlap=True, verbose=True)
np.save(out_dir / "coulomb_mat_overlap.npy", coulomb_overlap)

temperature_spec = 10  # Kelvin

lever_arm_matrix = np.load(out_dir / "lever_arm_matrix.npy")

gate_labels = ["BL", "P1", "P2", "BR"]

many_body_solver_params.energies = energies
many_body_solver_params.overlap = True
many_body_solver_params.coulomb_mat = coulomb_overlap
many_body_solver_params.alpha = lever_arm_matrix

jc = Junction(
    many_body_solver_params=many_body_solver_params,
    temperature=temperature_spec,
    contact_labels=gate_labels,
)


def get_add_spectrum(junc, gate_labels, gate_biases, temperature, verbose=True):

    junc.set_biases(gate_labels, gate_biases, verbose=verbose)
    out = add_spectrum(junc, temperature=temperature)
    return out


# set the bounds for the sweep
gate_1_min = -0.2
gate_1_max = 0.5
gate_2_min = -0.2
gate_2_max = 0.5

gate_1_biases = np.linspace(gate_1_min, gate_1_max, 90)
gate_2_biases = np.linspace(gate_2_min, gate_2_max, 90)

add_spectrum_mat = np.zeros((len(gate_1_biases), len(gate_2_biases)))

for idx_1, gate_1_bias in enumerate(gate_1_biases):
    for idx_2, gate_2_bias in enumerate(gate_2_biases):
        add_spectrum_mat[idx_1, idx_2] = get_add_spectrum(
            jc,
            gate_labels,
            np.array([0, gate_1_bias, gate_2_bias, 0]),
            temperature_spec,
            verbose=False,
        )

        np.savetxt(out_dir / "addition_spectrum.txt", add_spectrum_mat)

add_spec_to_plot = np.flip(np.transpose(add_spectrum_mat), axis=0)

fig, axs = plt.subplots()
axs.set_xlabel("$V_{g1}$ (V)", fontsize=16)
axs.set_ylabel("$V_{g2}$ (V)", fontsize=16)

diff_conds_map = axs.imshow(
    add_spec_to_plot / np.max(add_spec_to_plot),
    cmap="jet",
    interpolation="bilinear",
    extent=[gate_1_min - 0.6, gate_1_max - 0.6, gate_2_min - 0.6, gate_2_max - 0.6],
    aspect="auto",
)

fig.colorbar(diff_conds_map, ax=axs, label="Response (arb. units)")
fig.tight_layout()
fig.savefig(out_dir / "charge_stability_diagram.png")
plt.show()
