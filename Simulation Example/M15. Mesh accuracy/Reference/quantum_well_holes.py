__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

from qtcad.device import constants as ct
from qtcad.device import materials as mt
from qtcad.device.mesh1d import Mesh
from qtcad.device import Device
from qtcad.device.schrodinger import Solver, SolverParams
import numpy as np
from matplotlib import pyplot as plt
import pathlib

# Paths
path = pathlib.Path(__file__).parent.resolve()
path_out = path / "output"
path_out.mkdir(exist_ok=True)
path_mesh = str(path / "meshes" / "quantum_well.msh")

# Load the mesh
scale = 1e-9
mesh = Mesh(scale, path_mesh)

# Create the device
# In this example we will consider a 6 valence band model
dvc = Device(mesh, conf_carriers="h", hole_kp_model="luttinger_kohn_foreman_6band")

# Set the material to be GaAs
mat = mt.GaAs
dvc.new_region("well", mat)

# Then create insulator boundaries
dvc.new_insulator("left_barrier")
dvc.new_insulator("right_barrier")

# Solve Schrodinger's equation
# Parameters
params = SolverParams()
params.num_states = 20
params.ti_directions = ["x", "y"]
params.tol = 1e-10
# Solver
s = Solver(dvc, solver_params=params)

# Energies at (0, 0)
s.solve()
print("Energies at (0,0)")
dvc.print_energies()


# Energies at (1, 1)/aB
s.solve(k=np.array([1, 1]) / ct.aB, mixing=True)
print("Energies at (1, 1)/aB")
dvc.print_energies()

# Band structure
K = np.linspace(0, 0.3 / ct.aB, 30)
basis_vec = np.array([1, 0])

# With mixing
band_structure = s.band_structure(K, basis_vec)

fig = plt.figure(figsize=(7, 4))
ax1 = fig.add_subplot(1, 1, 1)
ax1.plot(K * ct.aB, band_structure / ct.e, "-")
ax1.set_xlabel("$k_x$ ($a_B^{-1}$)")
ax1.set_ylabel("Energy (eV)")
ax1.set_ylim([0, 5])
ax1.set_title("Band structure including band mixing")
plt.show()
fig.savefig(path_out / "band_structure_with_mixing.png", dpi=150)
np.savetxt(
    path_out / "band_structure_with_mixing.csv",
    np.column_stack([K] + [band_structure[:, i] for i in range(band_structure.shape[1])]),
    header="K (1/m), E_band_0..N (J)",
)

# Without mixing
# Compute band structure using QTCAD
band_structure = s.band_structure(K, basis_vec, mixing=False)


# Analytic solutions
# Particle in a box energies
def particle_in_box_E(m, L, n):
    return n**2 * np.pi**2 * ct.hbar**2 / (2 * m * L**2)


g1, g2, g3, Delta = tuple(mat.hole_kp_params)
# Out of plane effective masses
mHH = ct.me / (g1 - 2 * g2)  # HH effective mass GaAs
mLH = ct.me / (g1 + 2 * g2)  # LH effective mass
mSO = ct.me / g1  # SO effective mass
# In plane effective masses
mHHip = ct.me / (g1 + g2)  # HH in-plane effective mass GaAs
mLHip = ct.me / (g1 - g2)  # LH in-plane effective mass GaAs
mSOip = ct.me / (g1)  # SO in-plane effective mass GaAs

# Analytic energies
PIB_HH = np.array([particle_in_box_E(mHH, 2e-9, n) for n in range(1, 5)])
free_HH = ct.hbar**2 / 2 / mHHip * K**2
E_HH = PIB_HH[np.newaxis, :] + free_HH[:, np.newaxis]

PIB_LH = np.array([particle_in_box_E(mLH, 2e-9, n) for n in range(1, 5)])
free_LH = ct.hbar**2 / 2 / mLHip * K**2
E_LH = PIB_LH[np.newaxis, :] + free_LH[:, np.newaxis]

PIB_SO = np.array([particle_in_box_E(mSO, 2e-9, n) for n in range(1, 5)])
free_SO = ct.hbar**2 / 2 / mSOip * K**2
E_SO = (PIB_SO[np.newaxis, :] + free_SO[:, np.newaxis]) + Delta

# Plot
fig = plt.figure(figsize=(7, 4))
ax1 = fig.add_subplot(1, 1, 1)
ax1.plot(K * ct.aB, band_structure / ct.e, "-")
ax1.plot(K * ct.aB, E_HH / ct.e, "--b")  # Plot analytic HH energies in blue
ax1.plot(K * ct.aB, E_LH / ct.e, "--r")  # Plot analytic LH energies in blue
ax1.plot(K * ct.aB, E_SO / ct.e, "--g")  # Plot analytic SO energies in blue
ax1.set_xlabel("$k_x$ ($a_B^{-1}$)")
ax1.set_ylabel("Energy (eV)")
ax1.set_ylim([0, 5])
ax1.set_title("Band structure in the absence of band mixing")
plt.show()
fig.savefig(path_out / "band_structure_no_mixing_vs_analytic.png", dpi=150)
np.savez(
    path_out / "band_structure_no_mixing_vs_analytic.npz",
    K=K,
    band_structure=band_structure,
    E_HH=E_HH,
    E_LH=E_LH,
    E_SO=E_SO,
    mHH=mHH,
    mLH=mLH,
    mSO=mSO,
    mHHip=mHHip,
    mLHip=mLHip,
    mSOip=mSOip,
)
