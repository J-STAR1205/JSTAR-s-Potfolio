__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import pathlib
import numpy as np
import matplotlib.pyplot as plt
from qtcad.device import io
from qtcad.device import constants as ct
from qtcad.device.mesh3d import SubMesh
from qtcad.device import SubDevice
from double_dot_fdsoi import get_double_dot_fdsoi, V_set
from chemical_potential import compute_chemical_potential

# Files -----------------------------------------------------------------------
script_dir = pathlib.Path(__file__).parent.resolve()

path_out = script_dir / "output"
path_local_out = path_out / "coulomb_peaks"
phi_in_file = path_out / "phi" / "phi.hdf5"

leverarm_file = str(path_local_out / "lever_arm.txt")
lever_arm_plot_file = str(path_local_out / "lever_arm_plot.svg")

# Setup ----------------------------------------------------------------------

# Create the device
d, set_region, qd_region = get_double_dot_fdsoi()
# load the potential
d.set_potential(io.load(phi_in_file, var_name="phi"))
# Create sub-device for the SET region
set_mesh = SubMesh(d.mesh, set_region)
set_device = SubDevice(d, set_mesh)
set_device.set_insulator_boundaries()

# V_set values
V_SET = np.array([V_set, V_set - 0.05])
# Number of particles
Nmin = 12
Nmax = 14
# Generate data
Mu = compute_chemical_potential(
    V_SET, Nmin, Nmax, d, set_device, out_path=path_local_out
)

# Coulomb peak positions ------------------------------------------------------

# Fit a line: y = alpha*x + b
lever_arm = np.zeros(Mu.shape[1])
intercept = np.zeros(Mu.shape[1])
for n in range(Mu.shape[1]):
    alpha, b = np.polyfit(V_SET, Mu[:, n] / ct.e, 1)
    lever_arm[n] = alpha
    intercept[n] = b
    print(f"n={n + Nmin + 1}: Lever Arm [eV/V] = {alpha}, Intercept [eV] = {b}")

# Compute Coulomb peak positions
# (x-intercept)
V_peaks = -intercept / lever_arm
print("Coulomb Peak Positions [V]:")
print(V_peaks)
# (Distance between peaks)
Delta_V_peaks = np.diff(V_peaks)
print("Distance between Coulomb Peaks [V]:")
print(Delta_V_peaks)

# Save
header = "Lever Arm [eV/V]    Intercept [eV]    Coulomb Peak Position [V]"
np.savetxt(
    leverarm_file,
    np.column_stack((lever_arm, intercept, V_peaks)),
    fmt="%.10f",
    header=header,
)

# Plot
plt.figure(figsize=(8, 5))
# (Add data to plot)
DeltaV = np.max(np.insert(V_SET, 0, np.max(V_peaks))) - np.min(
    np.insert(V_SET, 0, np.min(V_peaks))
)
V_line = np.linspace(min(V_SET) - DeltaV, max(V_SET) + DeltaV, 200)
plt.axhline(y=0, color="black", linewidth=3)
for j, alpha in enumerate(lever_arm):
    b = intercept[j]  # Intercept for given value of N
    mu_line = alpha * V_line + b
    plt.plot(V_line, mu_line, color="red", linewidth=2.0, linestyle="--")
    plt.axvline(V_peaks[j], color="black", linestyle="--", linewidth=2.0)
for j in range(Mu.shape[1]):
    plt.scatter(V_SET, Mu[:, j] / ct.e, marker="o", s=60, label=f"N = {j + Nmin + 1}")

# (Labels, legend, grid)
plt.xlabel("$V_{SET}$ [V]", fontsize=16)
plt.ylabel("$\\mu$ [eV]", fontsize=16)
plt.legend(fontsize=14)
plt.grid(True)
# (Save)
plt.tight_layout()
plt.savefig(lever_arm_plot_file, dpi=300)

# Compute addition energies ---------------------------------------------------
addition_energies = Mu[:, 1:] - Mu[:, :-1]
print("Addition Energies [eV]:")
print(addition_energies / ct.e)

# Charging energy
E_C = addition_energies[0]

# Capacitance between SET gate and the SET ----------------------------------
C_gate_set = -lever_arm / (E_C / ct.e) * ct.e * 1e18  # in attoFarads [aF]
print("Capacitance between SET gate and the SET [aF]:")
print(C_gate_set)
