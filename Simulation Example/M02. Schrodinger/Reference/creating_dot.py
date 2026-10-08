__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import numpy as np
from matplotlib import pyplot as plt
import pathlib
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device import constants as ct
from qtcad.device.analysis import linecut
from qtcad.device import io
from qtcad.device import materials as mt
from qtcad.device import Device, SubDevice
from qtcad.device.schrodinger import Solver, SolverParams

# Paths to mesh file and output files
script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes" / "lead_dot_lead.msh"
path_V = script_dir / "output" / "creating_dot_V.vtu"
path_psi0 = script_dir / "output" / "creating_dot_psi0.vtu"
path_psi1 = script_dir / "output" / "creating_dot_psi1.vtu"
path_vtm = script_dir / "output" / "creating_dot_data.vtm"

# Load and plot the mesh
scalingFactor = 1e-9
mesh = Mesh(scalingFactor, str(path_mesh))
mesh.show()

# Create the device from the mesh
d = Device(mesh, conf_carriers="e")

d.new_region("dot", mt.GaAs)
d.new_region("lead", mt.GaAs)
d.new_region("lead2", mt.GaAs)
d.new_region("barrier", mt.GaAs)
d.new_region("barrier2", mt.GaAs)

# Potential along y and z
# Position of the harmonic oscillator minimum
y0 = 0
z0 = 0
Ky = 2e-7  # Parabolic in y
Kz = 1e-4  # Parabolic in z


def V_func(x, y, z):
    Vy = Ky / 2.0 * (y - y0) ** 2
    Vz = Kz / 2.0 * (z - z0) ** 2
    return Vy + Vz


# Set the potential energy
d.set_V(V_func)

# Add the potential energy along x
d.add_to_V(12e-3 * ct.e, region="barrier")
d.add_to_V(12e-3 * ct.e, region="barrier2")
d.add_to_V(-10e-3 * ct.e, region="lead")
d.add_to_V(-10e-3 * ct.e, region="lead2")
d.add_to_V(-1e-2 * ct.e)

# Create submesh
dotmesh = SubMesh(d.mesh, ["dot", "barrier", "barrier2"])

# Create subdevice
dot = SubDevice(d, dotmesh)

# Solve the Schrodinger equation.
# # Configure the Schrödinger solver
params_schrod = SolverParams()
params_schrod.num_states = 3  # Number of energy levels to consider
params_schrod.tol = 1e-12  # Set the tolerance for convergence

# # Create solver
s = Solver(dot, solver_params=params_schrod)
s.solve()  # solve

dot.print_energies()

# Global nodes in dot region
xdot = dotmesh.glob_nodes[:, 0]
ydot = dotmesh.glob_nodes[:, 1]
zdot = dotmesh.glob_nodes[:, 2]

# Global nodes in full device
x = mesh.glob_nodes[:, 0]
y = mesh.glob_nodes[:, 1]
z = mesh.glob_nodes[:, 2]

# Linecut coordinates for wave function and potential energy
beginwf = (np.min(xdot), y0, z0)
endwf = (np.max(xdot), y0, z0)
beginV = (np.min(x), y0, z0)
endV = (np.max(x), y0, z0)

# Linecuts of ground state squared and potential energy
psiposx, psix0 = linecut(dotmesh, np.abs(dot.eigenfunctions[:, 0]) ** 2, beginwf, endwf)
psiposx += np.min(xdot)
Vposx, Vx = linecut(mesh, d.V, beginV, endV)
Vposx += np.min(x)

beginy = ((np.max(x) + np.min(x)) / 2, np.min(y), z0)
endy = ((np.max(x) + np.min(x)) / 2, np.max(y), z0)
psiposy, psiy0 = linecut(dotmesh, np.abs(dot.eigenfunctions[:, 0]) ** 2, beginy, endy)
psiposy += np.min(ydot)
Vposy, Vy = linecut(mesh, d.V, beginy, endy)
Vposy += np.min(y)

nm = 1e-9

# Plotting eigenstates along x
fig, ax1 = plt.subplots()
ax1.set_title("Linecut along x")
ax1.plot(psiposx / nm, psix0)
ax2 = ax1.twinx()
ax2.plot(Vposx / nm, Vx / ct.e, "--")
ax1.grid()
ax1.set_xlabel("x (nm)")
ax1.set_ylabel(r"$|\Psi(x, y_0, z_0)|^2 (1/m^3)$")
ax2.set_ylabel(r"$V (eV)$")

# Display
plt.tight_layout()
plt.show()

# Plotting eigenstates along y
fig, ax3 = plt.subplots()
ax3.set_title("Linecut along y")
ax3.plot(psiposy / nm, psiy0)
ax4 = ax3.twinx()
ax4.plot(Vposy / nm, Vy / ct.e, "--")
ax3.grid()
ax3.set_xlabel("y (nm)")
ax3.set_ylabel(r"$|\Psi(x_0 , y, z_0)|^2 (1/m^3)$")
ax4.set_ylabel(r"$V (eV)$")

# Display
plt.tight_layout()
plt.show()

# Saving in .vtu format
io.save(path_V, d.V, mesh)
io.save(path_psi0, dot.eigenfunctions[:, 0], dotmesh)
io.save(path_psi1, dot.eigenfunctions[:, 1], dotmesh)

# Saving in .vtm format
data_dict = {
    "ground state": dot.eigenfunctions[:, 0],
    "first excited state": dot.eigenfunctions[:, 1],
    "Potential": dot.phi,
    "Confinement Potential": dot.V,
}
io.save(path_vtm, data_dict, dotmesh)
