__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import pathlib
import numpy as np
from numpy.typing import NDArray
from qtcad.device import materials as mt
from qtcad.device import constants as ct
from qtcad.device.mesh3d import Mesh
from qtcad.device import Device, SubDevice
from qtcad.device.analysis import plot_slice

# Default values
V_set = 1.25  # Default SET potential (V)
V_qubit = 0.75  # Default qubit potential (V)
scaling = 1e-9  # scaling factor for the mesh
valley_splitting = 0.5e-3 * ct.e

# Visualization parameters
normal = (0.0, 0.0, 1.0)  # the plane of the slice is normal to the z axis
origin = (0.0, 0.0, -1.0 * scaling)  # the slice is inside the Si quantum well


def save_slice(
    device: Device | SubDevice, data: NDArray, file_name: str | pathlib.Path, label: str
) -> None:
    """Helper function to save a slice of the device data.

    Args:
        device: The Device object containing the data.
        data: The data array to be sliced and saved.
        file_name: The filename to save the slice plot.
        label: The label for the colorbar axis.
    """
    plot_slice(
        device.mesh,
        data,
        normal=normal,
        origin=origin,
        cb_axis_label=label,
        path=file_name,
        show_figure=False,
    )


def state_probability_density(device: Device | SubDevice, state: int) -> NDArray:
    """Return the scalar probability density associated with an eigenstate.

    For multiband solutions, such as the explicit valley-resolved electron
    states used in this FD-SOI example, the final axis indexes band, spin,
    and/or valley components. Summing
    :math:`\\sum_\\nu |\\psi_{i,\\nu}(\\mathbf{r})|^2` over that axis recovers
    the scalar density used for visualization and charge-density
    calculations.

    Args:
        device: Device containing previously computed eigenfunctions.
        state: Single-particle state index for which to compute the
            probability density.

    Returns:
        Scalar probability density evaluated on the device mesh.
    """
    density = np.abs(device.eigenfunctions[:, state]) ** 2
    if density.ndim > 1:
        density = np.sum(density, axis=-1)
    return density


# Function to create the FD-SOI device
def get_double_dot_fdsoi(
    mesh_file_name: str = "refined_dqdfdsoi.msh",
    V_set: float = V_set,
    V_qubit: float = V_qubit,
    valley_splitting: float = valley_splitting,
) -> tuple[Device, list[str], list[str]]:
    """Create a Device object for the FD-SOI structure.

    Args:
        mesh_file_name: Filename of the mesh file.
        V_set: Potential applied to the SET gate.
        V_qubit: Potential applied to the qubit gate.
        valley_splitting: Valley splitting applied to the device.

    Returns:
        Device: The constructed Device object.
        list[str]: List of region names associated with the SET.
        list[str]: List of region names associated with the qubit quantum-dot.
    """
    # Create the mesh
    script_dir = pathlib.Path(__file__).parent.resolve()
    path_mesh = script_dir / "meshes"
    mesh_file = path_mesh / mesh_file_name
    mesh = Mesh(scaling, mesh_file)

    # Define the device object
    dvc = Device(mesh, conf_carriers="e")
    dvc.set_temperature(0.1)

    # Create the regions
    dvc.new_region("oxide", mt.SiO2)
    dvc.new_region("oxide.QD1", mt.SiO2)
    dvc.new_region("oxide.QD2", mt.SiO2)
    dvc.new_region("channel", mt.Si)
    dvc.new_region("channel.QD1", mt.Si)
    dvc.new_region("channel.QD2", mt.Si)
    dvc.new_region("source", mt.Si, ndoping=1e20 * 1e6)
    dvc.new_region("drain", mt.Si, ndoping=1e20 * 1e6)

    # Set up boundary conditions
    Ew = mt.Si.Eg / 2 + mt.Si.chi  # Midgap
    dvc.new_gate_bnd("barrier_gate_1_bnd", 0.5, Ew)
    dvc.new_gate_bnd("plunger_gate_1_bnd", V_set, Ew)
    dvc.new_gate_bnd("barrier_gate_2_bnd", 0.5, Ew)
    dvc.new_gate_bnd("plunger_gate_2_bnd", V_qubit, Ew)
    dvc.new_gate_bnd("barrier_gate_3_bnd", 0.5, Ew)
    dvc.new_ohmic_bnd("source_bnd")
    dvc.new_ohmic_bnd("drain_bnd")
    dvc.new_frozen_bnd("back_gate_bnd", -0.5, mt.Si, 1e15 * 1e6, "n", 46 * 1e-3 * ct.e)

    # Create the double quantum dot region
    SET_region_list = ["oxide.QD1", "channel.QD1"]
    qubit_region_list = ["oxide.QD2", "channel.QD2"]

    # Add valley splitting
    dvc.set_valley_splitting(valley_splitting)

    return dvc, SET_region_list, qubit_region_list
