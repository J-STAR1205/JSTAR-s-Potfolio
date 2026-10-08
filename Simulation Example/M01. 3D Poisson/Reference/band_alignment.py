__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
from qtcad.device import constants as ct
from qtcad.device.mesh1d import Mesh
from qtcad.device import Device
from qtcad.device import materials as mt
from qtcad.device import analysis as an
from qtcad.device.schrodinger import Solver

# Path to mesh file
path = Path(__file__).parent.resolve()
path = path / "meshes" / "band_alignment.msh"

# Load the mesh
mesh = Mesh(1e-9, path)

# ----------------------------------------------------------------------
# GaAs/AlGaAs — Anderson’s rule, RESCU alignment, and Schrödinger solver
# ----------------------------------------------------------------------

# Create the device and add materials
dvc = Device(mesh, conf_carriers="e")

dvc.new_region("left_barrier", mt.AlGaAs)
dvc.new_region("well", mt.GaAs)
dvc.new_region("right_barrier", mt.AlGaAs)

# Show the band diagram obtained from Anderson's rule
an.plot_bands(dvc, title="Anderson's rule")

# View the total confinement potential before calling set_V_from_phi
an.plot(
    mesh,
    dvc.get_Vconf() / ct.e,
    ylabel="$V_\\mathrm{conf}$ (eV)",
    title="Before setting V from phi",
)

# View the total confinement potential after calling set_V_from_phi
dvc.set_V_from_phi()
an.plot(
    mesh,
    dvc.get_Vconf() / ct.e,
    ylabel="$V_\\mathrm{conf}$ (eV)",
    title="After setting V from phi",
)

# Set an external potential to shift the barriers by +0.5 eV
dvc.set_Vext(0.5 * ct.e, "left_barrier")
dvc.set_Vext(0.5 * ct.e, "right_barrier")
an.plot(
    mesh,
    dvc.get_Vconf() / ct.e,
    ylabel="$V_\\mathrm{conf}$ (eV)",
    title="Shifted barrier heights",
)

# Unset the external potential
dvc.set_Vext(0.0, "left_barrier")
dvc.set_Vext(0.0, "right_barrier")

# Use RESCU-fitted alignment parameters; choose GaAs as the reference layer
dvc.align_bands(mt.GaAs)
an.plot_bands(dvc, title="Band alignment from RESCU simulations")

# Calculate the electron envelope functions
dvc.set_V_from_phi()
slv = Solver(dvc)
slv.solve()
dvc.print_energies()

# Plot the first few levels
x = mesh.glob_nodes[:, 0]
sort_indices = np.argsort(x)  # Sort the nodes along x
fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel("$x$ (nm)")
ax.set_ylabel("$\\psi(x)$")
ax.plot(
    x[sort_indices] / 1e-9,
    dvc.eigenfunctions[sort_indices, 0],
    "-k",
    label="Ground state",
)
ax.plot(
    x[sort_indices] / 1e-9,
    dvc.eigenfunctions[sort_indices, 1],
    "--b",
    label="1st excited state",
)
ax.plot(
    x[sort_indices] / 1e-9,
    dvc.eigenfunctions[sort_indices, 2],
    ":r",
    label="2nd excited state",
)
ax.legend()
plt.show()

# Shift the bands in the well downwards by 0.5 eV
dvc.add_to_ref_potential(-0.5, region="well")
an.plot_bands(
    dvc, title="Band alignment from RESCU simulations with -0.5 eV shift in well"
)

# Update the confinement potential and solve Schrödinger again
dvc.set_V_from_phi()
slv = Solver(dvc)
slv.solve()
dvc.print_energies()

# Plot the first few levels again
x = mesh.glob_nodes[:, 0]
sort_indices = np.argsort(x)  # Sort the nodes along x
fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel("$x$ (nm)")
ax.set_ylabel("$\\psi(x)$")
ax.plot(
    x[sort_indices] / 1e-9,
    dvc.eigenfunctions[sort_indices, 0],
    "-k",
    label="Ground state",
)
ax.plot(
    x[sort_indices] / 1e-9,
    dvc.eigenfunctions[sort_indices, 1],
    "--b",
    label="1st excited state",
)
ax.plot(
    x[sort_indices] / 1e-9,
    dvc.eigenfunctions[sort_indices, 2],
    ":r",
    label="2nd excited state",
)
ax.legend()
plt.show()

# ----------------------------------------------
# Band alignment in SiGe/Si/SiGe heterostructure
# ----------------------------------------------

# materials for SiGe/Si/SiGe heterostructure
mt.SiGe_DFT.set_alloy_composition(0.30)  # relaxed Si0.70Ge0.30 barriers
mt.Si_strained_on_SiGe.set_alloy_composition(0.30)

# Create the device
dvc = Device(mesh, conf_carriers="e")
dvc.new_region("left_barrier", mt.SiGe_DFT)
dvc.new_region("well", mt.Si_strained_on_SiGe)
dvc.new_region("right_barrier", mt.SiGe_DFT)

# Align bands using the relaxed SiGe as the reference layer
dvc.align_bands(mt.SiGe_DFT)

an.plot_bands(dvc, title="Si/SiGe band alignment (RESCU)")

# Conduction-band and valence-band offsets
cb = dvc.cond_band_edge() / ct.e
cbo = np.ptp(np.unique(np.round(cb, 6)))
vb = dvc.vlnce_band_edge() / ct.e
vbo = np.ptp(np.unique(np.round(vb, 6)))
print(f"SiGe/Si/SiGe heterostructure: VBO = {vbo:.6f} eV | CBO = {cbo:.6f} eV")

dvc.set_V_from_phi()
an.plot(
    mesh,
    dvc.get_Vconf() / ct.e,
    ylabel="$V_{\\mathrm{conf}}$ (eV)",
    title="Confinement potential (electron well)",
)

# ----------------------------------------------
# Band alignment in SiGe/Ge/SiGe heterostructure
# ----------------------------------------------

# materials for SiGe/Ge/SiGe heterostructure
mt.SiGe_DFT.set_alloy_composition(0.75)  # relaxed Si0.25Ge0.75 barriers
mt.Ge_strained_on_SiGe.set_alloy_composition(0.75)

# Create the device
dvc = Device(mesh, conf_carriers="h", hole_kp_model="luttinger_kohn_foreman")
dvc.new_region("left_barrier", mt.SiGe_DFT)
dvc.new_region("well", mt.Ge_strained_on_SiGe)
dvc.new_region("right_barrier", mt.SiGe_DFT)
dvc.align_bands(mt.SiGe_DFT)

an.plot_bands(dvc, title="Ge/SiGe band alignment (RESCU)")

vb = dvc.vlnce_band_edge() / ct.e
vbo = np.ptp(np.unique(np.round(vb, 6)))
cb = dvc.cond_band_edge() / ct.e
cbo = np.ptp(np.unique(np.round(cb, 6)))
print(f"SiGe/Ge/SiGe heterostructure: VBO = {vbo:.6f} eV | CBO = {cbo:.6f} eV")

dvc.set_V_from_phi()
an.plot(
    mesh,
    dvc.get_Vconf() / ct.e,
    ylabel="$V_{\\mathrm{conf}}$ (eV)",
    title="Confinement potential (hole well)",
)

# -----------------------------------------------------------------------
# Offset sweep – varying Ge content in the SiGe substrates (SiGe/Si/SiGe)
# -----------------------------------------------------------------------
xs = np.linspace(0.0, 0.8, 17)
print("-" * 31)
print("     x     VBO (eV)    CBO (eV)")
print("-" * 31)

for xval in xs:
    # Composition-dependent materials
    mt.SiGe_DFT.set_alloy_composition(xval)
    mt.Si_strained_on_SiGe.set_alloy_composition(xval)

    # Device & alignment (electron well)
    dvc = Device(mesh, conf_carriers="e")
    dvc.new_region("left_barrier", mt.SiGe_DFT)
    dvc.new_region("well", mt.Si_strained_on_SiGe)
    dvc.new_region("right_barrier", mt.SiGe_DFT)
    dvc.align_bands(mt.SiGe_DFT)

    cb = dvc.cond_band_edge() / ct.e
    vb = dvc.vlnce_band_edge() / ct.e
    cbo = np.ptp(np.unique(np.round(cb, 6)))
    vbo = np.ptp(np.unique(np.round(vb, 6)))

    print(f"  {xval:0.3f}   {vbo:0.6f}    {cbo:0.6f}")
