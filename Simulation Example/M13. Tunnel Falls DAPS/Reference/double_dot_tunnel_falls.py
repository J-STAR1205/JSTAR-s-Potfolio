__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

from pathlib import Path
import numpy as np
from qtcad.device import Device, SubDevice
from qtcad.device import materials as mt
from qtcad.device.analysis import plot_slice
from qtcad.device.mesh3d import Mesh

# Default values --------------------------------------------------------------
nm = 1e-9
Ge_fraction = 0.30  # Ge fraction in the SiGe layers.

# Visualization parameters
buffer_thickness = 50.0
qw_thickness = 4.6
normal = (0.0, 0.0, 1.0)
origin = (0.0, 0.0, (buffer_thickness + 0.5 * qw_thickness) * nm)

# Default gate voltages
V_B4 = 1.0
V_P5 = 2.5
V_B5 = 1.2
V_P6 = 2.5
V_B6 = 1.0
V_SG = -0.5
V_CS = -0.5


def save_slice(
    device: Device | SubDevice,
    data: np.ndarray,
    file_name: str | Path,
    label: str,
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


def state_probability_density(device: Device | SubDevice, state: int) -> np.ndarray:
    """Return the scalar probability density associated with an eigenstate.

    Args:
        device: Device containing previously computed eigenfunctions.
        state: Single-particle state index for which to compute the density.

    Returns:
        Scalar probability density evaluated on the device mesh.
    """
    density = np.abs(device.eigenfunctions[:, state]) ** 2
    if density.ndim > 1:
        density = np.sum(density, axis=-1)
    return density


def get_double_dot_tunnel_falls(
    mesh_file_name: str = "tunnel_falls_double_dot.msh",
    V_B4: float = V_B4,
    V_P5: float = V_P5,
    V_B5: float = V_B5,
    V_P6: float = V_P6,
    V_B6: float = V_B6,
    V_SG: float = V_SG,
    V_CS: float = V_CS,
) -> tuple[Device, list[str], list[str]]:
    """Create a Device object for the simplified Tunnel Falls double-dot mesh.

    Args:
        mesh_file_name: Name of the mesh file stored in the local ``meshes``
            directory.
        V_B4: Potential applied to gate ``B4``.
        V_P5: Potential applied to gate ``P5``.
        V_B5: Potential applied to gate ``B5``.
        V_P6: Potential applied to gate ``P6``.
        V_B6: Potential applied to gate ``B6``.
        V_SG: Potential applied to gate ``SG``.
        V_CS: Potential applied to gate ``CS``.

    Returns:
        Device, left-dot region labels, and right-dot region labels.
    """

    script_dir = Path(__file__).parent.resolve()
    # Load the mesh.
    mesh_file = script_dir / "meshes" / mesh_file_name
    mesh = Mesh(nm, mesh_file)

    # Create the device.
    device = Device(mesh, conf_carriers="e")

    # Generate materials for quantum-well and barrier regions.
    sige_barrier = mt.SiGe_DFT.copy()
    sige_barrier.set_alloy_composition(Ge_fraction)

    strained_well = mt.Si_strained_on_SiGe.copy()
    strained_well.set_alloy_composition(Ge_fraction)

    # Create the regions.
    device.new_region("relaxed_buffer", sige_barrier)
    device.new_region("quantum_well", strained_well)
    device.new_region("upper_barrier", sige_barrier)
    device.new_region("si_cap", mt.Si)
    device.new_region("gate_oxide_sio2", mt.SiO2)
    device.new_region("gate_oxide_hfo2", mt.HfO2)
    device.new_region("screening_ild", mt.SiO2)

    device.new_region("relaxed_buffer.dot_left", sige_barrier)
    device.new_region("quantum_well.dot_left", strained_well)
    device.new_region("upper_barrier.dot_left", sige_barrier)
    device.new_region("relaxed_buffer.dot_right", sige_barrier)
    device.new_region("quantum_well.dot_right", strained_well)
    device.new_region("upper_barrier.dot_right", sige_barrier)

    # Set up boundary conditions
    gate_work_function = strained_well.Eg / 2 + strained_well.chi
    gate_biases = {
        "B4": V_B4,
        "P5": V_P5,
        "B5": V_B5,
        "P6": V_P6,
        "B6": V_B6,
        "SG": V_SG,
        "CS": V_CS,
    }
    for gate_name, gate_bias in gate_biases.items():
        device.new_gate_bnd(gate_name, gate_bias, gate_work_function)

    # Create the double quantum dot region.
    left_dot_region = [
        "relaxed_buffer.dot_left",
        "quantum_well.dot_left",
        "upper_barrier.dot_left",
    ]
    right_dot_region = [
        "relaxed_buffer.dot_right",
        "quantum_well.dot_right",
        "upper_barrier.dot_right",
    ]

    return device, left_dot_region, right_dot_region
