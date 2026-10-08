__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import numpy as np
from matplotlib import pyplot as plt
from qtcad.device import constants as ct
import qutip
from qtcad.qubit import spectra, dynamics, noise

# Copied here are the Hamiltonian and drive for the projected 2 level system
# analyzed in the example `MOS_EDSR.py`. The Hamiltonian and drive are stored in the
# variables h0 and u on line 145 of `MOS_EDSR.py`.
# System Hamiltonian
h0 = np.array([[1.1128812e-23, 0.00000000e00], [0.00000000e00, 0.00000000e00]])
H0 = h0 / ct.hbar  # Rewrite in units of hbar = 1
# Drive
u = np.array(
    [
        [-2.03067270e-29, 1.00561612e-27 - 1.10285789e-37j],
        [1.00561612e-27 + 1.10285789e-37j, 0.0 + 0.0j],
    ]
)
delta_V = u / ct.hbar  # Rewrite in units of hbar = 1

omega = np.absolute(H0[0, 0] - H0[1, 1])  # qubit frequency
omega_Rabi = np.absolute(delta_V[0, 1])  # Rabi frequency on resonance
T_Rabi = 2 * np.pi / omega_Rabi  # Rabi period

T = 20 * T_Rabi  # Total time of the simulations
times = np.linspace(0, T, 1000)  # Times at which dynamics are simulated.

# Initialize objects that will handle dynamics and noise.
dyn = dynamics.Dynamics()
N = noise.Noise()

# The Hamiltonian H = H0 + delta_V cos(omega*t) in a frame rotating
# at frequency omega.
H = dyn.H_RF(H0, delta_V, omega)
# Solve for the dynamics without any noise.
psi0 = qutip.basis(2, 0)  # initial state
result0 = qutip.mesolve(H, psi0, times)

# Plot expectation values of sigma_z and sigma_y as a function of time.
fig, axes = plt.subplots(2)
fig.set_size_inches((8, 8))
# <sigma_z (t)>
axes[0].plot(result0.times, qutip.expect(qutip.sigmaz(), result0.states))
# <sigma_y (t)>
axes[1].plot(result0.times, qutip.expect(qutip.sigmay(), result0.states))
# Labels
axes[0].set_ylabel(f"$\left< \sigma_z (t) \\right>$", fontsize=20)
axes[1].set_ylabel(f"$\left< \sigma_y (t) \\right>$", fontsize=20)
axes[0].set_xlabel(f"$t \, (s)$", fontsize=20)
axes[1].set_xlabel(f"$t \, (s)$", fontsize=20)
plt.show()

# Include charge noise which is assumed to affect only the gate which
# generates the drive delta_V.
# Lorentzian spectrum
S0 = 0.01  # Total power
wc = omega_Rabi  # Cutoff frequency
w0 = 0  # Central frequency
# T0 should be much larger than the largest time scale of the problem under consideration.
# Used to properly define the Fourier series of the time series related to the noise.
T0 = 50 * T_Rabi  # Must take T0 such that 2*pi/T0 << wc to resolve the spectrum
omega_max = 5 * wc
wk = np.arange(0, omega_max, 2 * np.pi / T0) + w0
spectrum = spectra.lorentz(wk, S0, wc, w0=w0)

# Plot spectrum
fig, axes = plt.subplots(1)
fig.set_size_inches((10, 6))
axes.plot(wk, spectrum)
axes.set_ylabel(f"$S_\\beta(\omega) \, (s)$", fontsize=20)
axes.set_xlabel("$\omega \, (s^{-1})$", fontsize=20)
axes.set_title("Lorentzian spectral density", fontsize=30)
plt.show()

# Plot time series for the charge noise affecting the gate
# responsible for delta_V.
beta = N.gen_process(times, spectrum, wk, T0)
fig, axes = plt.subplots(1)
fig.set_size_inches((8, 5))
axes.plot(times, beta)
axes.set_ylabel(f"$\\beta(t)$", fontsize=20)
axes.set_xlabel("$t \, (s)$", fontsize=20)
plt.show()

# Simulate dynamic including noise.
num_runs1 = 1  # number of runs to average over
sigz_1run = N.dynamics(
    H0,
    delta_V,
    omega,
    spectrum,
    T0,
    psi0,
    times,
    qutip.sigmaz(),
    vec_omega=wk,
    num_runs=num_runs1,
)
num_runs2 = 100  # number of runs to average over
sigz_multiruns = N.dynamics(
    H0,
    delta_V,
    omega,
    spectrum,
    T0,
    psi0,
    times,
    qutip.sigmaz(),
    vec_omega=wk,
    num_runs=num_runs2,
)

# Plot dynamics
fig, axes = plt.subplots(2)
axes[0].plot(
    result0.times, qutip.expect(qutip.sigmaz(), result0.states), label="Without noise"
)
axes[0].plot(result0.times, sigz_1run, label="With noise")
axes[1].plot(
    result0.times, qutip.expect(qutip.sigmaz(), result0.states), label="Without noise"
)
axes[1].plot(result0.times, sigz_multiruns, label="With noise")
# Format Plots
axes[0].set_ylabel(f"$\left< \sigma_z (t) \\right>$", fontsize=20)
axes[1].set_ylabel(f"$\left< \sigma_z (t) \\right>$", fontsize=20)
axes[0].set_xlabel(f"$t \, (s)$", fontsize=20)
axes[1].set_xlabel(f"$t \, (s)$", fontsize=20)
axes[0].set_title(f"num_runs = {num_runs1}")
axes[1].set_title(f"num_runs = {num_runs2}")
axes[0].legend(loc="center left", fontsize=16)
axes[1].legend(loc="center left", fontsize=16)
fig.tight_layout()

plt.show()
