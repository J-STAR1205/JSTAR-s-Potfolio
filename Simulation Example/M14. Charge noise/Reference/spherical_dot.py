__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from qtcad.device import constants as ct
from qtcad.device.mesh3d import Mesh
from qtcad.device.materials import Si
from qtcad.atoms import Atoms, SubAtoms
from qtcad.atoms.unit_cell import UnitCellZincblende
from qtcad.atoms.keating import Solver as KeatingSolver
from qtcad.atoms.keating import SolverParams as KeatingSolverParams
from qtcad.atoms.schrodinger import Solver as SchrodingerSolver
from qtcad.atoms.schrodinger import SolverParams as SchrodingerSolverParams
from qtcad.atoms.analysis import save_vtu, analyze_dot

# Geometry of the spherical Si quantum dot in a Ge matrix
diameter = 4e-9
box_size = diameter * 2
x = np.array([-box_size / 2, box_size / 2])
y = x.copy()
z = x.copy()

# Unit cell of bulk Ge
lattice_constant = 5.647e-10
unit_cell = UnitCellZincblende(
    lattice_constant=lattice_constant, atom_species_ws=np.array(["Ge", "Ge"])
)

# Construct atomic structure from mesh
script_dir = Path(__file__).parent.resolve()
(script_dir / "output").mkdir(exist_ok=True)
path_mesh = script_dir / "meshes" / "spherical_dot.msh"
mesh = Mesh(scaling_factor=diameter, filename=str(path_mesh))
atoms_from_mesh = Atoms(mesh=mesh, unit_cell=unit_cell)
atoms_from_mesh.new_region(atom_species=np.array(["Si"] * 8), tag="ball")
path_vtu = script_dir / "output" / "spherical_dot_from_mesh.vtu"
save_vtu(atoms=atoms_from_mesh, out_dict={}, path=path_vtu)

# Axes along which periodic boundary conditions are applied for the atomic
# structure relaxation
pbc_directions = ["x", "y", "z"]

# Create atomic structure, defining the Si spherical dot via the `is_inside`
# argument of the `new_region` method
atoms = Atoms(x=x, y=y, z=z, unit_cell=unit_cell, pbc_directions=pbc_directions)


def is_inside_dot(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    return r <= diameter / 2


atoms.new_region(atom_species=np.array(["Si"] * 8), is_inside=is_inside_dot)

# Relax atomic structure with periodic boundary conditions
keating_solver_params = KeatingSolverParams()
keating_solver_params.verbose = True
keating_solver = KeatingSolver(atoms=atoms, solver_params=keating_solver_params)
keating_solver.solve()


# Impose electric potential corresponding to an electric field of 10 V/μm
# pointing in the -z direction
def phi(x, y, z):
    return 1e7 * z


atoms.set_potential(phi=phi)

# Define geometry of box where the Schrödinger equation is solved
box_size_tb = diameter * 1.5
x_tb = np.array([-box_size_tb / 2, box_size_tb / 2])
y_tb = x_tb.copy()
z_tb = x_tb.copy()

# Solve the Schrödinger equation without periodic boundary conditions
subatoms = SubAtoms(parent=atoms, x=x_tb, y=y_tb, z=z_tb, pbc_directions=[])
schrodinger_solver_params = SchrodingerSolverParams()
schrodinger_solver_params.verbose = True
schrodinger_solver_params.num_states = 5
schrodinger_solver = SchrodingerSolver(
    atoms=subatoms, solver_params=schrodinger_solver_params
)
energy_target = 0.7 * ct.e
schrodinger_solver.solve(energy_target=energy_target)
subatoms.print_energies()
analyze_dot(atoms=subatoms, psi=subatoms.eigenfunctions[0], verbose=True)
vtu_name = "spherical_dot.vtu"
vtu_path = script_dir / "output" / vtu_name
out_dict = {
    "Potential": subatoms.get_potential(),
    "Ground state": subatoms.eigenfunctions[0].get_prob_on_atoms(),
}
save_vtu(atoms=subatoms, out_dict=out_dict, path=vtu_path)

# Define and plot potential describing a positive charge defect at the origin
eps = Si.eps


def phi_defect(x, y, z):
    R = 1e-10
    r = np.sqrt(x**2 + y**2 + z**2)
    if r < R:
        return ct.e / (4 * np.pi * eps * R)
    else:
        return ct.e / (4 * np.pi * eps * r)


r_num = 100
r_vec = np.linspace(0, box_size_tb / 2, r_num)
phi_vec = np.array([phi_defect(r, 0, 0) for r in r_vec])
plt.figure()
plt.plot(r_vec * 1e9, phi_vec)
plt.xlabel("Distance from center of QD (nm)")
plt.ylabel("Electric potential (V)")
plt.title("Electric potential due to positive charge defect in spherical QD")
plt.grid()
plt.show()
plt.savefig(str(script_dir / "output" / "defect_potential_profile.png"), dpi=150)
np.savetxt(
    script_dir / "output" / "defect_potential_profile.csv",
    np.column_stack([r_vec, phi_vec]),
    header="r (m), phi_defect (V)",
)
atoms.add_to_potential(phi=phi_defect)

# Solve the Schrödinger equation with a positive charge defect at the center of
# the dot
subatoms_defect = SubAtoms(parent=atoms, x=x_tb, y=y_tb, z=z_tb, pbc_directions=[])
schrodinger_solver_defect = SchrodingerSolver(
    atoms=subatoms_defect, solver_params=schrodinger_solver_params
)
schrodinger_solver_defect.solve(energy_target=energy_target)
subatoms_defect.print_energies()
analyze_dot(atoms=subatoms_defect, psi=subatoms_defect.eigenfunctions[0], verbose=True)
vtu_name = "spherical_dot_defect.vtu"
vtu_path = script_dir / "output" / vtu_name
out_dict = {
    "Potential": subatoms_defect.get_potential(),
    "Ground state": subatoms_defect.eigenfunctions[0].get_prob_on_atoms(),
}
save_vtu(atoms=subatoms_defect, out_dict=out_dict, path=vtu_path)

# Define and plot functions specifying Si and Ge concentration profiles
diffusion_radius = diameter / 10


def Si_concentration(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    if r <= diameter / 2 - diffusion_radius:
        return 1.0
    elif r >= diameter / 2 + diffusion_radius:
        return 0.0
    else:
        a = -1 / (2 * diffusion_radius)
        b = (diameter / 2 + diffusion_radius) / (2 * diffusion_radius)
        return a * r + b


def Ge_concentration(x, y, z):
    return 1.0 - Si_concentration(x, y, z)


Si_vec = np.array([Si_concentration(r, 0, 0) for r in r_vec])
Ge_vec = np.array([Ge_concentration(r, 0, 0) for r in r_vec])
plt.figure()
plt.plot(r_vec * 1e9, Si_vec, label="Si")
plt.plot(r_vec * 1e9, Ge_vec, label="Ge")
plt.xlabel("Distance from center of QD (nm)")
plt.ylabel("Concentration of chemical species")
plt.title("Si and Ge concentration profiles in spherical QD")
plt.legend()
plt.grid()
plt.show()
plt.savefig(str(script_dir / "output" / "concentration_profiles.png"), dpi=150)
np.savetxt(
    script_dir / "output" / "concentration_profiles.csv",
    np.column_stack([r_vec, Si_vec, Ge_vec]),
    header="r (m), Si_concentration, Ge_concentration",
)

# Seed random number generator for reproducibility of the outputs of this
# tutorial; not required for simulations, in practice
rng_seed = 0

# Include gradients of Si and Ge concentration around in the dot-matrix
# interface
atoms_gradient = Atoms(
    x=x, y=y, z=z, unit_cell=unit_cell, pbc_directions=pbc_directions, rng_seed=rng_seed
)
atoms_gradient.new_region(
    atom_species=np.array([["Si", Si_concentration], ["Ge", Ge_concentration]])
)

# Relax atomic structure and solve the Schrödinger equation
keating_solver_gradient = KeatingSolver(
    atoms=atoms_gradient, solver_params=keating_solver_params
)
keating_solver_gradient.solve()
atoms_gradient.set_potential(phi=phi)
subatoms_gradient = SubAtoms(
    parent=atoms_gradient, x=x_tb, y=y_tb, z=z_tb, pbc_directions=[]
)
schrodinger_solver_gradient = SchrodingerSolver(
    atoms=subatoms_gradient, solver_params=schrodinger_solver_params
)
schrodinger_solver_gradient.solve(energy_target=energy_target)
subatoms_gradient.print_energies()
analyze_dot(
    atoms=subatoms_gradient, psi=subatoms_gradient.eigenfunctions[0], verbose=True
)
vtu_name = "spherical_dot_gradient.vtu"
vtu_path = script_dir / "output" / vtu_name
out_dict = {
    "Potential": subatoms_gradient.get_potential(),
    "Ground state": subatoms_gradient.eigenfunctions[0].get_prob_on_atoms(),
}
save_vtu(atoms=subatoms_gradient, out_dict=out_dict, path=vtu_path)
