__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import pathlib
import numpy as np
import matplotlib

# Use a non-interactive backend: the script only ever calls fig.savefig()
# (never plt.show()), so no display/event loop is needed. Forcing "Agg"
# here rules out any possibility of a GUI backend (Tk/Qt) stalling when
# run non-interactively / from an automated shell.
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import qtcad.device.constants as ct
from qtcad.device.mesh3d import Mesh
from qtcad.device import Device
from qtcad.device import io
from qtcad.atoms.atoms import Atoms
from qtcad.atoms.rough_surface import RoughSurface
from qtcad.atoms.schrodinger import Solver, SolverParams
from qtcad.atoms.analysis import save_vtu
from qtcad.atoms.g_tensor import Solver as GTensorSolver

# Setup
script_dir = pathlib.Path(__file__).parent.resolve()
out_dir = script_dir / "output"

## Load the mesh
path_mesh = script_dir / "meshes" / "refined_dqdfdsoi.msh"
scaling = 1e-9
mesh = Mesh(scaling, path_mesh)

## Create device and set potential in device to result of FEM simulation
d = Device(mesh, conf_carriers="e")
path_phi_hdf5 = out_dir / "detuned_potential.hdf5"
d.set_potential(io.load(path_phi_hdf5))

# Dimensions - Box below the left plunger gate
# NOTE: shrunk from the tutorial's original box (x:[-200,200]A, y:[-250,-100]A,
# z:[-80,10]A -> 269,934 atoms) to roughly 1/3.5 the atom count. The original
# size caused the tight-binding Schrodinger solver's shift-invert eigensolver
# to repeatedly fail to find eigenstates near energy_target and retry with a
# random perturbation (QTCAD's built-in robustness mechanism), with each
# attempt taking ~3-4 hours on this machine -- multiple retries made the run
# impractically long (7.7+ CPU-hours with only 2 of 10 allowed attempts done).
# This smaller box is a pragmatic reduction to get a working reference result;
# it no longer matches the tutorial's exact atom count/geometry extent.
A = 1e-10  # Angstroms
x = np.array([-120.0, 120.0]) * A
y = np.array([-220.0, -130.0]) * A
z = np.array([-60.0, 10.0]) * A

# Construct the atomic structure of the heterostructure
atoms = Atoms(x=x, y=y, z=z, verbose=True)

# Define the material stack of the heterostructure
## Create rough surface describing interface between Si and SiO2 layers
rng_seed = 1000
interface = RoughSurface(hurst=0.3, mean=0.0, rms=2.0 * A, rng_seed=rng_seed)
## Create the heterostructure layer
SiO2 = np.array([["SiO2", 1.0]])
atoms.new_layer(z_top=10.0 * A, z_bot=interface, atom_species=SiO2)

atoms.set_potential_from_device(d)
atoms.set_soc(True)
atoms.set_Bfield(np.array([0.0, 0.0, 1e-6]))
atoms.set_Zeeman(True)
atoms.set_tb_params(("Klimeck_cb", "Kim"))


# Define strain function
def strain_func(x, y, z):
    alpha_SiO2 = 0.55e-6
    alpha_Si = 2.6e-6
    delta_T = 1000.0 - 0.0  # K
    epsilon_parallel = (alpha_SiO2 - alpha_Si) * delta_T
    lambda_ = 2.0 * 1e-9  # 2 nm
    if z > 0:
        epsilon = epsilon_parallel
    else:
        epsilon = -epsilon_parallel * np.exp(z / lambda_)

    return np.array([[epsilon, 0.0, 0.0], [0.0, epsilon, 0.0], [0.0, 0.0, 0.0]])


# Set the strain
atoms.strain_atomic_structure(strain_func)

# Solve the Schrödinger equation

## Parametrize the Schrödinger solver
solver_params = SolverParams()
solver_params.num_states = 10
solver_params.verbose = True
# NOTE: memory_mode="robust" (the default) gives the best numerical
# stability, which this SOC/g-tensor problem seems to need (it has already
# failed and retried with a new random tight-binding perturbation several
# times across earlier attempts). But the default memory_threshold=0
# disables disk offloading entirely, so a retry that happens to need many
# more eigensolver iterations than usual can grow its working set without
# bound -- one run crashed system free memory from ~10 GB to ~270 MB on a
# single retry attempt that never finished. memory_mode="low" avoids that
# but is much slower even when memory is not actually tight (one run on
# this same (shrunk) atom box took 5.8+ CPU-hours in "low" mode vs ~10 min
# per attempt in "robust" mode). Setting memory_threshold=None keeps
# "robust" mode's speed/stability for the common case, but caps it at 75%
# of system RAM (~12 GiB here) as a safety net so a bad retry spills to
# disk instead of exhausting memory outright.
solver_params.memory_threshold = None
solver = Solver(atoms, solver_params=solver_params)

# Solve the Schrödinger equation
energy_target = 0.04 * ct.e
solver.solve(energy_target)
atoms.print_energies()

## Save modulus-square of projection of states on atoms to a .vtu file
psi = atoms.eigenfunctions
out_dict = {f"state {i}": psi[i].get_prob_on_atoms() for i in range(len(psi))}
path_vtu = out_dir / "FD-SOI_psi.vtu"
save_vtu(atoms=atoms, out_dict=out_dict, path=path_vtu, verbose=True)

# Compute the g-tensor
states = (psi[1], psi[0])
g_tensor_solver = GTensorSolver(atoms, states)
g = g_tensor_solver.solve()

print("g-tensor")
print(f"g = {g}")

## Zeeman splitting for in-plane magnetic field
Theta = np.linspace(0, np.pi, 41)  # angles in radians
Delta = []
for theta in Theta:
    B0 = 1  # units: T
    # Vary direction of in-plane field
    B = B0 * np.array([np.cos(theta), np.sin(theta), 0])
    Delta.append(g_tensor_solver.get_Zeeman_splitting(B))

Delta = np.array(Delta) / (ct.h * 1e6)  # Convert to MHz / T
Delta0 = np.mean(Delta)
print("Average spin-splitting: " + str(Delta0) + " MHz / T")

## Plotting the spin splitting
fig, ax = plt.subplots()
ax.plot(Theta / np.pi, (Delta - Delta0))
ax.set_xlabel(r"$\theta$ [$\pi$]")
ax.set_ylabel(r"$\Delta - \Delta_0$ [MHz / T]")
ax.grid()
fig.savefig(out_dir / "in_plane_g.png", dpi=300)

## Save spin splitting
with open(out_dir / "in_plane_g.txt", "a") as f:
    np.savetxt(f, Theta[np.newaxis, :] / np.pi, fmt="%.6f")
    np.savetxt(f, Delta[np.newaxis, :], fmt="%.6f")
