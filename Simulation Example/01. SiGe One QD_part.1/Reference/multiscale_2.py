__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import pathlib
import numpy as np
from qtcad.device import constants as ct
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device import Device, SubDevice
from qtcad.device import io
from qtcad.atoms.unit_cell import UnitCellZincblende
from qtcad.atoms import Atoms, SubAtoms
from qtcad.atoms.rough_surface import RoughSurface
from qtcad.atoms.keating import Solver as KeatingSolver
from qtcad.atoms.keating import SolverParams as KeatingSolverParams
from qtcad.atoms.schrodinger import Solver as SchrodingerSolver
from qtcad.atoms.schrodinger import SolverParams as SchrodingerSolverParams
from qtcad.atoms.analysis import (
    get_valley_phase,
    analyze_dot,
    save_vtu,
    get_operator_element,
    get_operator_from_fem,
)

# Load the mesh
script_dir = pathlib.Path(__file__).parent.resolve()
path_refined_mesh = str(script_dir / "meshes" / "multiscale.msh")
scaling = 1e-9
mesh = Mesh(scaling, path_refined_mesh)

# Create device and set potential in device to result of FEM simulation
d = Device(mesh, conf_carriers="e")
path_phi_hdf5 = str(script_dir / "output" / "multiscale_phi.hdf5")
d.set_potential(io.load(path_phi_hdf5, "phi"))

# Create unit cells of Si and Ge
unit_cell_Si = UnitCellZincblende(
    lattice_constant=5.431 * 1e-10, atom_species_ws=np.array(["Si", "Si"])
)
unit_cell_Ge = UnitCellZincblende(
    lattice_constant=5.658 * 1e-10, atom_species_ws=np.array(["Ge", "Ge"])
)

# Plot band structures of Si and Ge
path_bs_Si = str(script_dir / "output" / "multiscale_bs_Si.png")
unit_cell_Si.plot_bandstructure(title="Si", path=path_bs_Si, show_figure=False)
path_bs_Ge = str(script_dir / "output" / "multiscale_bs_Ge.png")
unit_cell_Ge.plot_bandstructure(title="Ge", path=path_bs_Ge, show_figure=False)

# Create atomic structure for heterostructure
x = np.array([-10.0, 10.0]) * 1e-9
y = np.array([-10.0, 10.0]) * 1e-9
z = np.array([-50.0, -5.0]) * 1e-9
atoms = Atoms(x=x, y=y, z=z, rng_seed=104)

# Create rough surface describing top and bottom surfaces of the well
well_top = RoughSurface(hurst=0.3, mean=-35.0 * 1e-9, rms=5.0 * 1e-10, rng_seed=0)
well_bot = RoughSurface(hurst=0.3, mean=-38.0 * 1e-9, rms=5.0 * 1e-10, rng_seed=1)

# Plot rough surface
x_plot = np.linspace(x[0], x[1], 100)
y_plot = np.linspace(y[0], y[1], 100)
path_rough_2d = str(script_dir / "output" / "multiscale_rough_2d.png")
well_top.plot(x=x_plot, y=y_plot, type="2D", path=path_rough_2d, show_figure=False)
path_rough_3d = str(script_dir / "output" / "multiscale_rough_3d.png")
well_top.plot(x=x_plot, y=y_plot, type="3D", path=path_rough_3d, show_figure=False)

# Define the material stack of the heterostructure
SiGe_alloy = np.array([["Si", 0.7], ["Ge", 0.3]])
atoms.new_layer(z_top=-5.0 * 1e-9, z_bot=well_top, atom_species=SiGe_alloy)  # cap
atoms.new_layer(z_top=well_bot, z_bot=z[0], atom_species=SiGe_alloy)  # buffer

# Relax atomic coordinates
keating_solver_params = KeatingSolverParams()
keating_solver_params.verbose = True
keating_solver = KeatingSolver(atoms=atoms, solver_params=keating_solver_params)
keating_solver.solve()
path_xyz = str(script_dir / "output" / "multiscale.xyz")
atoms.save_xyz(path=path_xyz)

# Contruct restricted atomic structure for QD
x_QD = np.array([-8.0, 8.0]) * 1e-9
y_QD = np.array([-8.0, 8.0]) * 1e-9
z_QD = np.array([-40.0, -33.0]) * 1e-9
atoms_QD = SubAtoms(parent=atoms, x=x_QD, y=y_QD, z=z_QD)

# Create Schrödinger equation solver within TB formalism
atoms_QD.set_potential_from_device(d=d)
schrodinger_solver_params = SchrodingerSolverParams()
schrodinger_solver_params.verbose = True
schrodinger_solver_params.num_states = 10
schrodinger_solver = SchrodingerSolver(
    atoms=atoms_QD, solver_params=schrodinger_solver_params
)

# Solve the Schrödinger equation
energy_target = -0.1 * ct.e
schrodinger_solver.solve(energy_target=energy_target)
atoms_QD.print_energies()

# Compute valley splitting
valley_splitting = atoms_QD.energies[1] - atoms_QD.energies[0]
print("Valley splitting: " + str(valley_splitting / ct.e) + " eV")

# Store valley states in variables nu_1 and nu_2
nu_1 = atoms_QD.eigenfunctions[0]
nu_2 = atoms_QD.eigenfunctions[1]

# Compute valley phase
valley_phase = get_valley_phase(atoms=atoms_QD, psi_0=nu_1, psi_1=nu_2)
print("Valley phase: " + str(valley_phase) + " radians")

# Compute modulus-square of projections of ground state on the
# spds* orbitals and on the chemical species in the system
# Also compute the geometric properties of the ground-state wavefunction
print("Ground state (nu_1)")
nu_1.get_prob_on_orbitals(verbose=True)
nu_1.get_prob_on_chemical_species(atoms=atoms_QD, verbose=True)
analyze_dot(atoms=atoms_QD, psi=nu_1, verbose=True)

# Repeat these steps for the first excited state
print("First excited state (nu_2)")
nu_2.get_prob_on_orbitals(verbose=True)
nu_2.get_prob_on_chemical_species(atoms=atoms_QD, verbose=True)
analyze_dot(atoms=atoms_QD, psi=nu_2, verbose=True)

# Save modulus-square of projection of states on atoms to a .vtu file
prob_on_atoms_0 = nu_1.get_prob_on_atoms()
prob_on_atoms_1 = nu_2.get_prob_on_atoms()
prob_on_atoms_2 = atoms_QD.eigenfunctions[2].get_prob_on_atoms()
path_vtu = str(script_dir / "output" / "multiscale_TB_wf.vtu")
out_dict = {
    "Probability density (ground state)": prob_on_atoms_0,
    "Probability density (1st excited state)": prob_on_atoms_1,
    "Probability density (2nd excited state)": prob_on_atoms_2,
}
save_vtu(atoms=atoms_QD, out_dict=out_dict, path=path_vtu)

# Compute the matrix element of the spin-orbit coupling Hamiltonian between
# the valley states with opposite spins
H_soc = schrodinger_solver.get_hamiltonian_soc()
nu_2_up = nu_2.add_spin_to_orbital_state(spin=True, in_place=False)
nu_1_down = nu_1.add_spin_to_orbital_state(spin=False, in_place=False)
C = get_operator_element(psi_l=nu_2_up, op=H_soc, psi_r=nu_1_down)
print("Matrix element of SOC Hamiltonian: " + str(C / ct.e) + " eV")

# Compute the matrix element of the derivative of the potential with respect to
# barrier gate voltage between the valley states
qd = ["cap_qd", "well_qd", "buffer_qd"]  # quantum dot regions
submesh = SubMesh(mesh, qd)
subd = SubDevice(d, submesh)
path_der_hdf5 = str(script_dir / "output" / "multiscale_der.hdf5")
der = io.load(path_der_hdf5, "der")
Dfg = get_operator_from_fem(
    atoms=atoms_QD, tb_orbs=nu_1.tb_orbs, tb_spin=nu_1.tb_spin, d=subd, quantity=der
)
D = get_operator_element(psi_l=nu_1, op=Dfg, psi_r=nu_2)
print(
    "Matrix element of derivative of potential with respect to barrier "
    + "gate voltage: "
    + str(D / ct.e)
    + " eV/V"
)

# Save parameters of quantum control problem to .txt file
params = np.array([valley_splitting, C, D])
path_params = str(script_dir / "output" / "multiscale_params.txt")
np.savetxt(path_params, params)
