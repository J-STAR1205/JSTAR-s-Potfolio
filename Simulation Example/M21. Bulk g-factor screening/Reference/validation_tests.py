__copyright__ = "Copyright 2026"

# M21 Bulk g-factor screening -- validation tests (spec section 17).
# Uses a trivial confined-electron box (NOT the real FD-SOI device yet) to
# verify the Device-level (non-atoms) Zeeman/g_star API behaves as expected
# before connecting it to the real device and running any sweep.
#
# This script intentionally does NOT import qtcad.atoms anywhere.

import pathlib
import numpy as np
from qtcad.device import constants as ct
from qtcad.device import materials as mt
from qtcad.device.mesh3d import Mesh
from qtcad.device import Device
from qtcad.device.schrodinger import Solver, SolverParams

script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes" / "box_order2.msh"

nm = 1e-9
mesh = Mesh(nm, str(path_mesh))

mu_B = 9.2740100783e-24  # Bohr magneton (J/T)
G_BULK = 1.998  # close to the free-electron value, used only for this test


def make_device(B_vec, g_star, zeeman=True):
    """Fresh device each time: simple parabolic confinement + constant B."""
    d = Device(mesh, conf_carriers="e")
    d.new_region("domain", mt.Si)
    d.new_insulator("bnd")
    x = mesh.glob_nodes[:, 0]
    y = mesh.glob_nodes[:, 1]
    z = mesh.glob_nodes[:, 2]
    K = 5e-3  # parabolic confinement strength (arbitrary, just needs bound states)
    V = 0.5 * K * ct.e * (x**2 + y**2 + z**2) / nm**2
    d.set_V(V)
    d.set_g_star(g_star)
    d.set_Zeeman(zeeman)
    if B_vec is not None:
        d.set_Bfield(np.array(B_vec))
    return d


def solve_lowest(d, num_states=4):
    params = SolverParams()
    params.num_states = num_states
    params.verbose = False
    s = Solver(d, solver_params=params)
    s.solve()
    return d.energies  # 1D array, ascending


print("=" * 80)
print("TEST A: B = 0 -> spin degeneracy should be preserved")
print("=" * 80)
dA = make_device(B_vec=[0.0, 0.0, 0.0], g_star=G_BULK)
EA = solve_lowest(dA)
print("Energies (eV):", EA / ct.e)
pair_gap_A = EA[1] - EA[0]
print(f"Gap between state 0 and 1: {pair_gap_A/ct.e*1e6:.6f} ueV "
      f"({'degenerate (OK)' if pair_gap_A < 1e-6*ct.e else 'NOT degenerate -- check setup'})")

print()
print("=" * 80)
print("TEST B: B = 0.1 T along x -> Zeeman splitting should appear")
print("=" * 80)
dB = make_device(B_vec=[0.1, 0.0, 0.0], g_star=G_BULK)
EB = solve_lowest(dB)
print("Energies (eV):", EB / ct.e)
dE_B = EB[1] - EB[0]
dE_B_analytic = G_BULK * mu_B * 0.1
print(f"QTCAD splitting   : {dE_B/ct.e*1e6:.6f} ueV")
print(f"Analytic g*mu_B*B : {dE_B_analytic/ct.e*1e6:.6f} ueV")
print(f"Relative error    : {(dE_B-dE_B_analytic)/dE_B_analytic*100:+.3f} %")

print()
print("=" * 80)
print("TEST C: B = 0.5 T along x -> splitting should scale ~5x vs Test B")
print("=" * 80)
dC = make_device(B_vec=[0.5, 0.0, 0.0], g_star=G_BULK)
EC = solve_lowest(dC)
dE_C = EC[1] - EC[0]
print(f"dE(0.5T) / dE(0.1T) = {dE_C/dE_B:.4f} (expected close to 5.0)")

print()
print("=" * 80)
print("TEST D: B = 0.5 T along x / y / z -> isotropic g* should give same splitting")
print("=" * 80)
results_xyz = {}
for axis, Bvec in zip("xyz", ([0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5])):
    d_ = make_device(B_vec=Bvec, g_star=G_BULK)
    E_ = solve_lowest(d_)
    dE_ = E_[1] - E_[0]
    results_xyz[axis] = dE_
    print(f"B along {axis}: splitting = {dE_/ct.e*1e6:.6f} ueV, g_eff = {dE_/(mu_B*0.5):.5f}")

spread = (max(results_xyz.values()) - min(results_xyz.values())) / np.mean(list(results_xyz.values())) * 100
print(f"Spread across x/y/z: {spread:.4f} % (should be ~0 for isotropic g*)")

print()
print("=" * 80)
print("TEST E: changing g_bulk should scale the splitting proportionally")
print("=" * 80)
for g_test in (1.0, 1.998, 4.0):
    d_ = make_device(B_vec=[0.5, 0, 0], g_star=g_test)
    E_ = solve_lowest(d_)
    dE_ = E_[1] - E_[0]
    print(f"g_bulk = {g_test:.3f} -> splitting = {dE_/ct.e*1e6:.6f} ueV, "
          f"dE/g = {dE_/g_test/ct.e*1e6:.6f} ueV (should be constant across rows)")

print()
print("All validation tests completed.")
