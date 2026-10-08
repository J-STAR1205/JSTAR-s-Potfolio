__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import pathlib
import numpy as np
from copy import deepcopy
from matplotlib import pyplot as plt
from qtcad.device import constants as ct
from qtcad.device.mesh3d import Mesh
from qtcad.device import analysis
from qtcad.device import materials as mt
from qtcad.device import Device
from qtcad.device.poisson import Solver, SolverParams

# Define some device parameters
Vgs = 2.0  # Gate-source bias
Ew = mt.Si.chi + mt.Si.Eg / 2  # Metal work function
Ltot = 20e-9  # Device height in m
radius = 2.5e-9  # Device radius in m

# Background charge parameters
surf_charge_dnsty = -ct.e * 5e17  # Surface charge density in C/m^2
vol_charge = -5 * ct.e  # Total charge in background volume charge density in C
vol_charge_spread = 1e-9  # Standard deviation of Gaussian charge density in m

# Paths to mesh file and output files
script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes" / "nanowire_background_charges.msh"
path_geo = script_dir / "meshes" / "nanowire_background_charges.geo_unrolled"
path_out = script_dir / "output"
path_out.mkdir(exist_ok=True)

# Load the mesh
scaling = 1e-9
mesh = Mesh(scaling, str(path_mesh))

# Create device and set temperature to 10 mK
dvc = Device(mesh, conf_carriers="e")
dvc.set_temperature(10e-3)

# Create regions first. The last added region takes priority for
# nodes that are shared by multiple regions.
# Make sure each nodes is assigned to a region, otherwise default
# (silicon) parameters are used
dvc.new_region("channel_ox", mt.SiO2)
dvc.new_region("source_ox", mt.SiO2)
dvc.new_region("drain_ox", mt.SiO2)
dvc.new_region("channel", mt.Si)
dvc.new_region("source", mt.Si, pdoping=5e20 * 1e6, ndoping=0)
dvc.new_region("drain", mt.Si, pdoping=5e20 * 1e6, ndoping=0)

# Then create boundaries
dvc.new_ohmic_bnd("source_bnd")
dvc.new_ohmic_bnd("drain_bnd")
dvc.new_gate_bnd("gate_bnd", Vgs, Ew)

# Keep a copy of the device in which no background charge is used
dvc_no_charge = deepcopy(dvc)

# Add surface charge at interface between silicon channel and SiO2
dvc.set_surf_charge_density("interface", surf_charge_dnsty)

# Create function for the volume charge density
# Take a Gaussian-shaped charge density in the SiO2 in proximity to the Si
# channel
x0 = 0
y0 = 2.0e-9
z0 = 10e-9


def gaussian(x, y, z):
    norm_factor = vol_charge / ((2 * np.pi) ** (3 / 2.0) * vol_charge_spread**3)
    numerator = (x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2
    denominator = 2 * vol_charge_spread**2
    return norm_factor * np.exp(-numerator / denominator)


dvc.set_vol_charge_density(gaussian)

# Configure the non-linear Poisson solver
params_poisson = SolverParams()
params_poisson.tol = 1e-5  # Convergence threshold (tolerance) for the error
params_poisson.initial_ref_factor = 0.1
params_poisson.final_ref_factor = 0.75
params_poisson.min_nodes = 20000

# Create an adaptive-mesh non-linear Poisson solver for each device
slv = Solver(dvc, solver_params=params_poisson, geo_file=path_geo)
slv_no_charge = Solver(dvc_no_charge, solver_params=params_poisson, geo_file=path_geo)


# Self-consistent solution for each device
print("=" * 80)
print("SOLVING NON-LINEAR POISSON WITHOUT BACKGROUND CHARGES")
print("=" * 80)
slv_no_charge.solve()

print("=" * 80)
print("SOLVING NON-LINEAR POISSON WITH BACKGROUND CHARGES")
print("=" * 80)
slv.solve()

# Evaluate conduction band edge along linecuts to compare the two cases

## Linecut along symmetry axis of the nanowire
z, band_edge_z = analysis.linecut(
    dvc.mesh, dvc.cond_band_edge() / ct.e, (0, 0, 0), (0, 0, Ltot)
)
z_no_charge, band_edge_z_no_charge = analysis.linecut(
    dvc_no_charge.mesh, dvc_no_charge.cond_band_edge() / ct.e, (0, 0, 0), (0, 0, Ltot)
)

## Linecut along y axis
y, band_edge_y = analysis.linecut(
    dvc.mesh, dvc.cond_band_edge() / ct.e, (0, -radius, Ltot / 2), (0, radius, Ltot / 2)
)
y_no_charge, band_edge_y_no_charge = analysis.linecut(
    dvc_no_charge.mesh,
    dvc_no_charge.cond_band_edge() / ct.e,
    (0, -radius, Ltot / 2),
    (0, radius, Ltot / 2),
)

# Produce plots of the linecuts
fig = plt.figure(figsize=(8, 7))

ax1 = fig.add_subplot(2, 1, 1)
ax1.set_xlabel("$z$ (nm)")
ax1.set_ylabel("$E_C$ (eV)")
ax1.plot(z / 1e-9, band_edge_z, "-b", label="With background charge")
ax1.plot(
    z_no_charge / 1e-9, band_edge_z_no_charge, "-r", label="Without background charge"
)
ax1.legend()

ax2 = fig.add_subplot(2, 1, 2)
ax2.set_xlabel("$y$ (nm)")
ax2.set_ylabel("$E_C$ (eV)")
ax2.plot(y / 1e-9, band_edge_y, "-b", label="With background charge")
ax2.plot(
    y_no_charge / 1e-9, band_edge_y_no_charge, "-r", label="Without background charge"
)
ax2.legend()
plt.show()
fig.savefig(str(path_out / "background_charges_linecuts.png"), dpi=150)

# Save linecut data for downstream analysis.
# NOTE: the "with charge" and "no charge" devices are adaptively meshed
# independently (different node distributions along the linecut), so the two
# curves have different lengths and cannot be column-stacked; save separately.
np.savetxt(
    path_out / "background_charges_linecut_z_with_charge.csv",
    np.column_stack([z, band_edge_z]),
    header="z (m), E_C with charge (eV)",
)
np.savetxt(
    path_out / "background_charges_linecut_z_no_charge.csv",
    np.column_stack([z_no_charge, band_edge_z_no_charge]),
    header="z_no_charge (m), E_C no charge (eV)",
)
np.savetxt(
    path_out / "background_charges_linecut_y_with_charge.csv",
    np.column_stack([y, band_edge_y]),
    header="y (m), E_C with charge (eV)",
)
np.savetxt(
    path_out / "background_charges_linecut_y_no_charge.csv",
    np.column_stack([y_no_charge, band_edge_y_no_charge]),
    header="y_no_charge (m), E_C no charge (eV)",
)

# Plot charge density in device with background charges along slices.
# NOTE: analysis.plot_slices() launches a full VTK/pyvista off-screen render
# that took ~900 s in this environment (same pathological slowness as
# Builder.view() encountered in M12/M13) for negligible visual gain over a
# direct scatter of the mesh-node charge density near the x=0 plane. We
# reconstruct an equivalent cross-section with matplotlib instead, which is
# >100x faster and trivially reproducible in the analysis notebook.
# dvc.rho is stored per tetrahedral element with one value per corner node
# (shape (Ntetra, 4)), not per global node, so we use element centroids (mean
# of the 4 corner-node positions) and the per-element mean of the 4 nodal
# values as the plotting coordinates/colors, rather than dvc.mesh.glob_nodes
# directly.
nodes = dvc.mesh.glob_nodes
connectivity = np.asarray(dvc.mesh.connectivity)
centroids = nodes[connectivity].mean(axis=1)  # (Ntetra, 3)
rho_cm3 = (dvc.rho / ct.e / 1e6).mean(axis=1)  # (Ntetra,)
slice_mask = np.abs(centroids[:, 0]) < 0.3e-9  # near x=0 plane
np.savez(
    path_out / "background_charges_rho_slice.npz",
    y=centroids[slice_mask, 1],
    z=centroids[slice_mask, 2],
    rho_cm3=rho_cm3[slice_mask],
)
fig2, ax = plt.subplots(figsize=(6, 7))
sc = ax.scatter(
    centroids[slice_mask, 1] / 1e-9,
    centroids[slice_mask, 2] / 1e-9,
    c=rho_cm3[slice_mask],
    cmap="viridis",
    s=8,
)
fig2.colorbar(sc, ax=ax, label=r"$\rho/e$ (cm$^{-3}$)")
ax.set_xlabel("$y$ (nm)")
ax.set_ylabel("$z$ (nm)")
ax.set_title("Charge density near $x=0$ plane")
plt.show()
fig2.savefig(str(path_out / "background_charges_rho_slice.png"), dpi=150)
