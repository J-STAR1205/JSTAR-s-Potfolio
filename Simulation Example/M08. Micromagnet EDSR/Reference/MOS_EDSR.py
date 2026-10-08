__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import pathlib
import numpy as np
from matplotlib import pyplot as plt
from qtcad.device.mesh3d import Mesh
from qtcad.device import constants as ct
from qtcad.device import materials as mt
from qtcad.device import Device
from qtcad.device.schrodinger import Solver as ssolver
from qtcad.device.schrodinger import SolverParams as SchrodingerSolverParams
from qtcad.device.poisson import Solver as psolver
from qtcad.device.poisson import SolverParams as PoissonSolverParams
from qtcad.device import analysis
from qtcad.device.operators import Gate
from qtcad.qubit import dynamics
import qutip

# Mesh
# Path to mesh file
script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes/" / "MOS_EDSR_example.msh"

# Load the mesh
scaling = 1e-9
mesh = Mesh(scaling, str(path_mesh))

# Create device from mesh and set statistics.
# It is important to set statistics before creating boundaries since
# they determine the boundary values, Failing to do so will produce
# inconsistent results.
# Calculate electronic charge using wavefunctions
d = Device(mesh, conf_carriers="e")
# Analytic approximation to Fermi-Dirac statistics
d.statistics = "FD_approx"

# Create regions first. The last added region takes priority for
# nodes that are shared by multiple regions.
# Make sure each node is assigned to a region, otherwise default
# (silicon) parameters are used
# Barrier
d.new_region("barrier", mt.SiO2)
# Confined
d.new_region("confined", mt.Si, pdoping=0, ndoping=0)
# The rest
d.new_region("undoped", mt.Si, pdoping=0, ndoping=0)
d.new_region("doped", mt.Si, pdoping=5e18 * 1e6, ndoping=0)

# Create boundaries
d.new_ohmic_bnd("bottom")

phi_tg = 1.0  # Top-gate voltage
phi_sg = 0.0  # Side-gate voltage
Ew = mt.Si.chi + mt.Si.Eg / 2  # Metal work function
d.new_gate_bnd("top_gate_bnd", phi_tg, Ew)
d.new_gate_bnd("side_gate_bnd", phi_sg, Ew)

ps = psolver(d)
ps.solve()  # Solve the Poisson equation using default parameters

# Convert the electric potential from the Poisson solver into a potential
# for the Schrödinger solver
d.set_V_from_phi()


def Bfield(x, y, z):
    B0 = 0.6  # homogeneous magnetic field (T)
    b = 0.3 * 1e6  # magnetic-field gradient (T/micron)

    return np.array([B0, 0, b * x])


d.set_Bfield(Bfield)


# Configure the Schrödinger solver
params_schrod = SchrodingerSolverParams()
num_states = 6
params_schrod.num_states = num_states  # Number of energy levels to consider
params_schrod.tol = 1e-12
ss = ssolver(d, solver_params=params_schrod)

ss.solve()  # Solve the Schrödinger equation

states = np.sum(np.abs(d.eigenfunctions) ** 2, axis=2)
# Plot eigenstates
analysis.plot_slices(
    mesh, states[:, 0], x=0, y=0, z=30 * scaling, title="Ground-state density"
)
analysis.plot_slices(
    mesh, states[:, 1], x=0, y=0, z=30 * scaling, title="First-excited-state density"
)
analysis.plot_slices(
    mesh, states[:, 2], x=0, y=0, z=30 * scaling, title="Second-excited-state density"
)
analysis.plot_slices(
    mesh, states[:, 3], x=0, y=0, z=30 * scaling, title="Third-excited-state density"
)
analysis.plot_slices(
    mesh, states[:, 4], x=0, y=0, z=30 * scaling, title="Fourth-excited-state density"
)
analysis.plot_slices(
    mesh, states[:, 5], x=0, y=0, z=30 * scaling, title="Fifth-excited-state density"
)

# Energies
d.print_energies()

gate_biases = np.linspace(0, 2, num=6)
gate = "side_gate_bnd"
UU = np.zeros((len(gate_biases), num_states, num_states), dtype=complex)

# Parametrize Poisson solver (inside Gate operator)
gate_params = PoissonSolverParams()
gate_params.tol = 1e-3

for i, V in enumerate(gate_biases):
    # Apply V to gate and compute the potential energy operator matrix
    G = Gate(d, gate, V, params=gate_params)
    UU[i, :, :] = G.get_operator_matrix()

# The potential energy matrix is a linear function of the applied voltage
fig = plt.figure()

# Plot - 01 matrix element
ax1 = fig.add_subplot(111)
ax1.plot(gate_biases, np.abs(np.real(UU[:, 0, 1]) / ct.e), ".-", label="Real")
ax1.plot(gate_biases, np.abs(np.imag(UU[:, 0, 1]) / ct.e), ".-", label="Imaginary")
ax1.set_xlabel(r"$\varphi^\mathrm{bias}$", fontsize=20)
ax1.set_ylabel(r"$\delta V_{01}(\varphi^\mathrm{bias})$ (eV)", fontsize=20)
ax1.grid()
ax1.legend()
fig.tight_layout()
plt.show()

# Use the fact that the perturbing potential is a linear function of the
# applied voltage.
delta_V = UU[-1]

# System Hamiltonian in the eigenbasis is diagonal
E = d.energies
E = E - (E[0])  # Go into frame rotating at ground-state frequency

# Create Hamiltonian
H0 = np.diagflat(E)
print(f"evals (eV): {np.diagonal(H0) / ct.e}")

# Simple EDSR calculation. The calculation is performed by projecting onto
# the relevant 2 dimensional subspace. Here we are studying transitions
# from the ground state (j=0) to the first excited state (i=1). The optional
# parameter plot is set to true to visualize the results.
dyn = dynamics.Dynamics()
result, h0, u, omega0, omega_rabi = dyn.transition_2_levels(
    1, 0, H0, delta_V, plot=True, npts=4000
)

# [VALIDATION] [수정 — M08 Step 1, stale-data pipeline fix]
# 2-레벨 Hamiltonian h0, drive u, qubit/Rabi 각주파수(rad/s), 그리고 이 결과를 만든
# 파라미터(gate bias, B-field)를 명시적으로 저장한다. EDSR_noise.py는 이제 이 파일을
# 읽어서 쓰도록 수정되었다 — 과거에는 이 값들을 손으로 복사해 하드코딩했는데,
# MOS_EDSR.py의 파라미터가 바뀐 뒤에도 EDSR_noise.py가 갱신되지 않아 Rabi 주파수가
# 1.6899MHz(현재 실행) vs 1.5177MHz(과거 하드코딩)로 약 11% 불일치하는 문제가 있었음
# (M03_M06_M08_Physics_Review_Report.md 7절 참고). 단위 표기: omega0/omega_rabi는
# rad/s (각주파수), f_rabi_Hz는 Hz (omega_rabi / 2pi) — 두 표현을 혼용하지 않기 위해
# 둘 다 명시적으로 저장한다 ([NUMERICAL] 2pi 변환을 암묵적으로 하지 않음).
import json
import sys
import qtcad

output_dir = script_dir / "output"
output_dir.mkdir(exist_ok=True)
f_rabi_hz = omega_rabi / (2 * np.pi)
qubit_frequency_rad_s = omega0  # transition_2_levels()가 반환하는 0-1 전이각주파수
np.savez(
    output_dir / "mos_edsr_h0_u.npz",
    h0=h0,
    u=u,
    omega0=omega0,
    omega_rabi=omega_rabi,
    f_rabi_hz=f_rabi_hz,
    gate_biases=gate_biases,
    phi_tg=phi_tg,
    phi_sg=phi_sg,
    B0=0.6,
    b_gradient=0.3e6,
    num_states=num_states,
)
metadata = {
    "qtcad_version": getattr(qtcad, "__version__", "unknown (2.2.6 per CLAUDE.md)"),
    "python_version": sys.version,
    "source_script": "MOS_EDSR.py",
    "mesh_file": str(path_mesh),
    "phi_tg_V": phi_tg,
    "phi_sg_V": phi_sg,
    "B0_T": 0.6,
    "b_gradient_T_per_m": 0.3e6,
    "num_states": num_states,
    "transition_levels": [1, 0],
    "omega0_rad_s": float(omega0),
    "omega_rabi_rad_s": float(omega_rabi),
    "f_rabi_Hz": float(f_rabi_hz),
    "units_note": "omega0/omega_rabi in rad/s; f_rabi_Hz = omega_rabi/(2*pi) in Hz; "
                   "h0,u in J (joules).",
}
with open(output_dir / "mos_edsr_metadata.json", "w", encoding="utf-8") as f:
    json.dump(metadata, f, indent=2)
print(f"\n[저장] h0, u, omega0={omega0:.6e} rad/s, omega_rabi={omega_rabi:.6e} rad/s "
      f"-> {output_dir / 'mos_edsr_h0_u.npz'}")
print(f"[저장] Rabi frequency (Hz): {f_rabi_hz:.6e} -> {output_dir / 'mos_edsr_metadata.json'}")

# We simulate the dynamics of
# an electron initialized in the ground state subject to a modulating voltage
# applied to the side gate without projecting onto a two-dimensional subspace.
# This section of the code requires a basic understanding of the QuTip library
H0 = qutip.Qobj(H0 / ct.hbar)  # system Hamiltonian in units of hbar
# Reset the zero of energy (ground state has E!=0).
np.fill_diagonal(delta_V, delta_V.diagonal() - delta_V[0, 0])
delta_V = qutip.Qobj(
    delta_V / ct.hbar
)  # potential energy matrix amplitude in units of hbar
H = [
    H0,
    [delta_V, "cos(omega0*t)"],
]  # Full Hamiltonian - sinusoidal modulation of delta_V

args = {"omega0": omega0}

psi0 = qutip.basis(6, 0)  # initialize in the ground state

T_Rabi = 2 * np.pi / omega_rabi  # Rabi cycle in the two-level limit
# array of times for which the solver should store the state vector
tlist = np.linspace(0, 2 * T_Rabi, 500000)
print("Computing Rabi oscillations beyond two-dimensional subspace.")
print("This may take several minutes.")
result = qutip.mesolve(
    H,
    psi0,
    tlist,
    c_ops=[],
    e_ops=[],
    args=args,
    options={"progress_bar": "tqdm"},
)

# Use the Dynamics object to get projectors onto the 6 basis states.
Proj = dyn.projectors(6)
# Plot projections onto each state.
fig, axes = plt.subplots(2, 3)
for i in range(2):
    for j in range(3):
        axes[i, j].plot(tlist, qutip.expect(Proj[3 * i + j], result.states))
        axes[i, j].set_xlabel(r"$t$", fontsize=20)
        axes[i, j].set_ylabel(f"$P_{3 * i + j}$", fontsize=20)
fig.tight_layout()
plt.show()
