__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import pathlib
import numpy as np
from matplotlib import pyplot as plt
import qtcad.device.constants as ct

# Load the quantum control parameters from the TB simulation
script_dir = pathlib.Path(__file__).parent.resolve()
path_params = str(script_dir / "output" / "multiscale_params.txt")
params = np.loadtxt(path_params, dtype=complex)
V = np.real(params[0])  # valley splitting
C = params[1]  # matrix element of SOC Hamiltonian
D = np.real(params[2])  # matrix element of derivative of potential

# Related model parameters
E1 = -V / 2  # energy level of first valley state
E2 = V / 2  # energy level of second valley state
g = 2.0  # gyromagnetic factor
BA = V / (g * ct.muB)  # anticrossing magnetic field

# Sampled values of amplitude of magnetic field
Bmax = 1.5 * BA
Bnum = 1000
B = np.linspace(0.0, Bmax, Bnum)

# Parametrized energy levels of the quantum dot
E_QD_0 = E1 - 0.5 * g * ct.muB * B
E_QD_1 = 0.5 * (E1 + E2) - 0.5 * np.sqrt((V - g * ct.muB * B) ** 2 + 4 * np.abs(C) ** 2)
E_QD_2 = 0.5 * (E1 + E2) + 0.5 * np.sqrt((V - g * ct.muB * B) ** 2 + 4 * np.abs(C) ** 2)
E_QD_3 = E2 + 0.5 * g * ct.muB * B

# Plot the energy levels
plt.plot(B, 1e3 * E_QD_0 / ct.e, label=r"$\left|0\right>$")
plt.plot(B, 1e3 * E_QD_1 / ct.e, label=r"$\left|1\right>$")
plt.plot(B, 1e3 * E_QD_2 / ct.e, label=r"$\left|2\right>$")
plt.plot(B, 1e3 * E_QD_3 / ct.e, label=r"$\left|3\right>$")
plt.xlabel("Magnetic field (T)")
plt.ylabel("Energy (meV)")
plt.legend()
plt.tight_layout()
plt.show()

# Parametrized EDSR Rabi frequency
eps = V - g * ct.muB * B + np.sqrt((V - g * ct.muB * B) ** 2 + 4 * np.abs(C) ** 2)
alpha = 2 * np.conj(C) / np.sqrt(eps**2 + 4 * np.abs(C) ** 2)
beta = eps / np.sqrt(eps**2 + 4 * np.abs(C) ** 2)
epsp = V + g * ct.muB * B + np.sqrt((V + g * ct.muB * B) ** 2 + 4 * np.abs(C) ** 2)
alphap = -np.conj(2 * np.conj(C) / np.sqrt(epsp**2 + 4 * np.abs(C) ** 2))
betap = epsp / np.sqrt(epsp**2 + 4 * np.abs(C) ** 2)
VBG = 0.001  # amplitude of rf oscillations of barrier gate voltage
hf = VBG * np.abs(np.conj(alpha) * betap + alphap * beta) * np.abs(D)

# Plot the Rabi frequency
plt.plot(B, hf / ct.h / 1e6)
plt.xlabel("Magnetic field (T)")
plt.ylabel("Rabi frequency (MHz)")
plt.tight_layout()
plt.show()
