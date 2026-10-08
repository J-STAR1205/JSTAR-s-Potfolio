__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import numpy as np
import matplotlib.pyplot as plt
import pathlib
import os
from joblib import load
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device.device import Device, SubDevice
from qtcad.device.poisson import Solver, SolverParams
from qtcad.device.schrodinger import Solver as SchrodingerSolver
from qtcad.device.schrodinger import SolverParams as SchrodingerSolverParams
from qtcad.device import constants as ct
from helper.double_dot_fdsoi import get_double_dot_fdsoi

# Files -----------------------------------------------------------------------
script_dir = pathlib.Path(__file__).parent.resolve()

path_mesh = script_dir / "meshes"
path_out = script_dir / "output"

# NOTE: gmsh's XAO reader/writer fails on non-ASCII (Korean) paths with
# "Could not load XML file" (same bug hit in M12/M13). The project directory
# is under a Korean-named folder, so we stage the mesh/geometry files through
# an ASCII-only path in C:\temp and read/write there instead.
path_mesh_ascii = pathlib.Path(r"C:\temp\m14_dqdfdsoi")
path_mesh_ascii.mkdir(parents=True, exist_ok=True)

mesh_file = str(path_mesh / "dqdfdsoi.msh")
geo_file = str(path_mesh_ascii / "dqdfdsoi.xao")
ref_file = str(path_mesh_ascii / "refined_dqdfdsoi_oxide.msh")
donor_profile_file = str(path_out / "pc_density_profile_0.050_4.0.joblib")

E_file_G = str(path_out / "fdsoi_energies_oxide_charge_G.txt")
Delta_file_G = str(path_out / "DeltaE_oxide_G.png")
E_file_D = str(path_out / "fdsoi_energies_oxide_charge_D.txt")
Delta_file_D = str(path_out / "DeltaE_oxide_D.png")

# Setup -----------------------------------------------------------------------

# Load the mesh
scaling = 1e-9
mesh = Mesh(scaling, mesh_file)

# Define the gate bias parameters
back_gate_bias = -0.5
barrier_gate_1_bias = 0.5
plunger_gate_1_bias = 0.7
barrier_gate_2_bias = 0.5
plunger_gate_2_bias = 0.6
barrier_gate_3_bias = 0.5


# Poisson ---------------------------------------------------------------------
def func_poisson(pos: np.ndarray, Q: np.ndarray, r: float, prof="Gaussian") -> Device:
    """Apply the adaptive non-linear Poisson solver with point charges.
    Args:
        pos (np.ndarray): Position of the point charges
        Q (np.ndarray): Charge of the point charges
        r (float): Smearing radius of the point charges
        prof (str): Profile of the point charges
    Returns:
        device: Device over which the non-linear Poisson solver has been
            applied (point charges have been included).
    """

    # Define the device object from the function defined in the FD-SOI tutorial
    dvc = get_double_dot_fdsoi(
        mesh,
        back_gate_bias,
        barrier_gate_1_bias,
        plunger_gate_1_bias,
        barrier_gate_2_bias,
        plunger_gate_2_bias,
        barrier_gate_3_bias,
    )

    dvc.add_point_charges(pos, Q, r=r, refine=True, profile=prof)

    # Configure the Non-Linear Poisson solver
    solver_params = SolverParams()
    solver_params.tol = 1e-3
    solver_params.initial_ref_factor = 0.1
    solver_params.final_ref_factor = 0.75
    solver_params.min_nodes = 50000
    solver_params.max_nodes = 1e5
    solver_params.maxiter_adapt = 30
    solver_params.maxiter = 200
    solver_params.refined_region = [
        "oxide_dot",
        "gate_oxide_dot",
        "buried_oxide_dot",
        "channel_dot",
    ]
    solver_params.h_refined = 0.8
    solver_params.refined_mesh_filename = ref_file

    # Solve Poisson
    slv = Solver(dvc, solver_params=solver_params, geo_file=geo_file)
    slv.solve()

    return dvc


def func_schrodinger(d: Device) -> SubDevice:
    """Apply the Schrodinger solver to the double-dot region of the FD-SOI
    device.

    Args:
        d (Device): An FD-SOI device over which the non-linear Poisson
            solver has been applied.

    Returns:
        SubDevice: The subdevice object over which the Schödinger equation
            has been solved.
    """

    d.set_V_from_phi()

    # List of regions forming the double quantum dot region
    dot_region_list = ["oxide_dot", "gate_oxide_dot", "buried_oxide_dot", "channel_dot"]
    # Create a submesh including only the dot region
    submesh = SubMesh(d.mesh, dot_region_list)
    # Create a subdevice object
    subdevice = SubDevice(d, submesh)

    # Configure the Schrodinger solver
    solver_params = SchrodingerSolverParams()
    solver_params.num_states = 10
    solver_params.tol = 1e-6

    # Solve Schrodinger
    slv = SchrodingerSolver(subdevice, solver_params=solver_params)
    slv.solve()

    return subdevice


# Generate data ---------------------------------------------------------------


def generate_data(prof, save_file) -> None:
    """Compute the eigenergies of the FD-SOI as a function of the position
    of a point charge.

    Args:
        prof (str): Profile of the point charge.
        save_file (str): File to save the data.

    """

    # Loop over point-charge positions
    npts_dir = 3  # Number of points per direction
    x0, y0 = 0, -17.5  # Center of plunger gate 1
    xmin, ymin = -20, -25  # Minimum position

    x_var = np.linspace(xmin, x0, npts_dir)
    x_const = x0 * np.ones(npts_dir)
    y_var = np.linspace(ymin, y0, npts_dir)
    y_const = y0 * np.ones(npts_dir)

    x_values = np.append(x_var[0 : npts_dir - 1], x_const)
    y_values = np.append(y_const[0 : npts_dir - 1], y_var)

    for i, x in enumerate(x_values):
        y = y_values[i]

        pos = np.array([[x, y, 0]]) * scaling  # position of the point charge
        Q = np.array([1]) * ct.e  # charge of the point charge
        r = 1e-10  # radius of the point charge

        dvc = func_poisson(pos, Q, r, prof=prof)  # Solve Poisson
        subdevice = func_schrodinger(dvc)  # Solve Schrodinger

        subdevice.print_energies()  # Print energies

        # Save energies
        out = np.append(np.array([x, y]), subdevice.energies / ct.e)
        # header
        E_header = ""  # set empty header
        if not os.path.isfile(save_file):  # checks if the file exists
            # if it doesn't then add the header
            E_header = "x (nm), y (nm), energies (eV)"
        with open(save_file, "a") as file:
            np.savetxt(file, np.array(out)[np.newaxis, :], header=E_header)


# Visualize data --------------------------------------------------------------


def visualize_data(E_file, Delta_file) -> None:
    """Visualize the data generated by the generate_data function.

    Args:
        E_file (str): File containing the data.
        Delta_file (str): File to save the plot.
    """

    # Load the data
    data = np.loadtxt(E_file)
    x = data[:, 0]
    y = data[:, 1]
    gap = data[:, 3] - data[:, 2]

    # Plot labels
    label = "$\Delta E$ (eV)"
    title = "Orbital level spacing"
    xlabel = "$x$ (nm)"
    ylabel = "$y$ (nm)"

    # Create figure and axes
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot the data
    scatter = ax.scatter(x, y, c=gap, cmap="viridis")
    fig.colorbar(scatter, ax=ax, label=label)

    # Set labels
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True)

    # Display the plot
    plt.savefig(str(Delta_file))
    plt.show()


if __name__ == "__main__":
    # Generate data
    profile = "Gaussian"
    generate_data(profile, E_file_G)

    # Visualize data
    visualize_data(E_file_G, Delta_file_G)

    # The donor-profile branch below requires a precomputed donor charge
    # density (pc_density_profile_0.050_4.0.joblib) produced by the
    # donor_MVEMT.py tutorial (Device 13). Per CURRICULUM.md, Device 13 is
    # explicitly excluded from this curriculum ("제외한 튜토리얼"), so that
    # file was never generated and this branch is skipped. The Gaussian
    # point-charge profile above already exercises the full Poisson+
    # Schrodinger point-charge pipeline this module is meant to teach.
    if os.path.isfile(donor_profile_file):
        donor = load(donor_profile_file)
        # Generate data
        profile = donor
        generate_data(profile, E_file_D)

        # Visualize data
        visualize_data(E_file_D, Delta_file_D)
    else:
        print(
            "Skipping donor-profile branch: "
            f"{donor_profile_file} not found (Device 13 excluded, see CURRICULUM.md)."
        )
