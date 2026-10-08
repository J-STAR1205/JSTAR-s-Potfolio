__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from qtcad.device import Device, SubDevice
from qtcad.device import constants as ct
from qtcad.device.mesh3d import SubMesh
from qtcad.device.poisson_linear import Solver as LinearPoissonSolver
from qtcad.device.poisson_linear import SolverParams as LinearPoissonSolverParams
from qtcad.device.schrodinger import Solver as SchrodingerSolver
from qtcad.device.schrodinger import SolverParams as SchrodingerSolverParams
from double_dot_tunnel_falls import (
    get_double_dot_tunnel_falls,
    save_slice,
    state_probability_density,
)
from double_dot_tunnel_falls import V_B4, V_P5, V_B5, V_B6, V_SG, V_CS
from valley_kp_model import generate_correlated_two_valley_model

# File paths.
script_dir = Path(__file__).parent.resolve()
mesh_dir = script_dir / "meshes"
path_out = script_dir / "output"
path_local_out = path_out / "energy_vs_detuning"
path_local_out.mkdir(parents=True, exist_ok=True)

mesh_file = mesh_dir / "tunnel_falls_double_dot.msh"
# NOTE: gmsh fails to load .xao files for adaptive meshing on non-ASCII
# (Korean) Windows paths ("Could not load XML file" / similar). Use the
# ASCII-path copy written during the M13 Builder step instead.
geo_file = Path(r"C:\temp\m13_tunnel_falls\tunnel_falls_double_dot.xao")
refined_mesh_file = mesh_dir / "refined_tunnel_falls_double_dot.msh"

valley_map_file = path_local_out / "valley_splitting_map.png"
energy_spectrum_file = path_local_out / "energy_spectrum.txt"
energy_spectrum_plot_file = path_local_out / "energy_spectrum.png"
ground_state_density_files = [
    path_local_out / "ground_state_density_negative_detuning.png",
    path_local_out / "ground_state_density_positive_detuning.png",
]

# Detuning sweep.
V_P_center = V_P5
Delta_V = 0.01
num_voltages = 11
detuning_values = np.linspace(-Delta_V, Delta_V, num_voltages)

# Correlated valley-splitting parameters.
mean_valley_splitting = 200.0e-6 * ct.e
correlation_length_x = 19.2e-9
correlation_length_y = 1000e-9
map_nx = 200
map_ny = 15
random_seed = 8


def gate_biases_for_detuning(detuning: float) -> dict[str, float]:
    """Return the gate voltages used for one detuning point.

    Args:
        detuning: Symmetric P5-P6 detuning in volts. Positive detuning raises
            ``P5`` and lowers ``P6`` by the same amount.

    Returns:
        Dictionary of gate voltages accepted by the linear Poisson solver.
    """
    return {
        "V_B4": V_B4,
        "V_P5": V_P_center + 0.5 * detuning,
        "V_B5": V_B5,
        "V_P6": V_P_center - 0.5 * detuning,
        "V_B6": V_B6,
        "V_SG": V_SG,
        "V_CS": V_CS,
    }


def solve_linear_poisson(
    gate_biases: dict[str, float],
    *,
    mesh_file_name: str | Path = mesh_file,
    geo_file_name: str | Path | None = None,
    refined_mesh_file: str | Path | None = None,
) -> tuple[Device, list[str], list[str]]:
    """Solve the linear Poisson equation for one gate configuration.

    Args:
        gate_biases: Dictionary of gate voltages accepted by the Poisson solve.
        mesh_file_name: Mesh file used to create the device.
        geo_file_name: Geometry file used to activate adaptive meshing. When
            omitted, the Poisson equation is solved on the input mesh without
            adaptive meshing.
        refined_mesh_file: Output path for the refined mesh written by the
            adaptive solve.

    Returns:
        Device with solved potential, left-dot region labels, and right-dot
        region labels.
    """

    # Create the device.
    device, left_dot_region, right_dot_region = get_double_dot_tunnel_falls(
        mesh_file_name=str(mesh_file_name),
        **gate_biases,
    )

    # Set up linear Poisson solver.
    solver_params = LinearPoissonSolverParams()
    solver_params.tol = 1e-12
    # Adaptive meshing if a geometry file is passed.
    if geo_file_name is not None:
        solver_params.eta0 = 0.30
        solver_params.refined_region = left_dot_region + right_dot_region
        solver_params.h_refined = 1.0
        if refined_mesh_file is not None:
            solver_params.refined_mesh_filename = str(refined_mesh_file)
    solver = LinearPoissonSolver(
        d=device,
        solver_params=solver_params,
        geo_file=geo_file_name,
    )

    solver.solve()
    device.set_V_from_phi()

    return device, left_dot_region, right_dot_region


def solve_schrodinger(
    subdevice: Device | SubDevice,
    *,
    num_states: int = 8,
    tol: float = 1e-12,
) -> Device | SubDevice:
    """Solve the single-particle Schrödinger equation on a dot subdevice.

    Args:
        subdevice: Dot subdevice on which the solve is performed.
        num_states: Number of eigenstates requested from the solver.
        tol: Schrödinger solver tolerance.

    Returns:
        The solved subdevice.
    """
    # Schrödinger solver parameters.
    solver_params = SchrodingerSolverParams()
    solver_params.num_states = num_states
    solver_params.tol = tol

    # Force wavefunctions to 0 at the edges.
    subdevice.set_insulator_boundaries()

    solver = SchrodingerSolver(subdevice, solver_params=solver_params)
    solver.solve()
    return subdevice


def save_energy_spectrum_plot(detuning: np.ndarray, energy_data: np.ndarray) -> None:
    """Save the energy spectrum as a function of detuning.

    Args:
        detuning: Detuning values.
        energy_data: Double-dot eigenenergies for different voltage detuning values.
            Each row corresponds to an entry in detuning.
    """
    energies_mev = energy_data[:, :4] * 1e3
    energies_mev -= np.nanmin(energies_mev[:, 0])

    figure, axis = plt.subplots(figsize=(6.0, 4.5))
    for state in range(energies_mev.shape[1]):
        axis.plot(
            detuning * 1000,
            energies_mev[:, state],
            marker="o",
            linewidth=1.5,
            label=f"State {state}",
        )

    axis.set_xlabel("Detuning [mV]")
    axis.set_ylabel("$E - \\min (E_0)$ [meV]")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8, ncol=2)
    figure.tight_layout()
    figure.savefig(energy_spectrum_plot_file, dpi=200)
    plt.close(figure)


if __name__ == "__main__":
    # Generate the refined mesh.
    device, left_dot_region, right_dot_region = solve_linear_poisson(
        gate_biases_for_detuning(0.0),
        mesh_file_name=mesh_file,
        geo_file_name=geo_file,
        refined_mesh_file=refined_mesh_file,
    )

    # Get k·p model on the refined mesh.
    dot_region = left_dot_region + right_dot_region
    dot_mesh = SubMesh(device.mesh, dot_region)

    model = generate_correlated_two_valley_model(
        dot_mesh.glob_nodes[:, :2],
        mean_valley_splitting=mean_valley_splitting,
        correlation_length_x=correlation_length_x,
        correlation_length_y=correlation_length_y,
        map_nx=map_nx,
        map_ny=map_ny,
        random_seed=random_seed,
        valley_map_filename=valley_map_file,
    )

    energy_data = []

    for detuning_index, detuning in enumerate(detuning_values):
        gate_biases = gate_biases_for_detuning(detuning)

        print("\n" + "=" * 60)
        print(f"Detuning point {detuning_index + 1}/{detuning_values.size}")
        print(f"detuning = {detuning:.6f} V")
        print("=" * 60)

        # Solve the linear Poisson equation.
        device, _, _ = solve_linear_poisson(
            gate_biases,
            mesh_file_name=refined_mesh_file,
        )

        # Solve the Schrödinger equation.
        dot_device = SubDevice(device, dot_mesh)
        dot_device.set_electron_kp_model(model)
        solve_schrodinger(dot_device)
        dot_device.print_energies()

        # Save ground-state density slices at the detuning-sweep endpoints.
        if detuning_index in (0, detuning_values.size - 1):
            file_index = 0 if detuning_index == 0 else 1
            save_slice(
                dot_device,
                state_probability_density(dot_device, 0),
                ground_state_density_files[file_index],
                label=r"$|\psi_0|^2$ [1 / m$^3$]",
            )

        # Save energies.
        header = "Detuning [V]    Energies [eV]"
        with open(energy_spectrum_file, "a") as f:
            out = [detuning] + list(dot_device.energies / ct.e)
            np.savetxt(f, np.array(out)[np.newaxis, :], fmt="%.8f", header=header)

        energy_data.append(dot_device.energies / ct.e)

    # Plot spectrum.
    energy_data = np.asarray(energy_data)
    save_energy_spectrum_plot(detuning_values, energy_data)
