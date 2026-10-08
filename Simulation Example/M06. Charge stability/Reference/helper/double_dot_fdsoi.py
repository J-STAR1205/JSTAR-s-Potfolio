__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

from qtcad.device import constants as ct
from qtcad.device import materials as mt
from qtcad.device import Device


def get_double_dot_fdsoi(
    mesh,
    back_gate_bias,
    barrier_gate_1_bias,
    plunger_gate_1_bias,
    barrier_gate_2_bias,
    plunger_gate_2_bias,
    barrier_gate_3_bias,
):
    """Produce a Device object for a double quantum dot in a Fully-Depleted
    Silicon-On-Insulator (FD-SOI) transistor.

    Args:
       mesh (Mesh): The mesh over which the device is defined.
       back_gate_bias (float): The potential applied at the back gate.
       barrier_gate_1_bias (float): The potential applied at the first
          barrier gate.
       plunger_gate_1_bias (float): The potential applied at the first
          plunger gate.
       barrier_gate_2_bias (float): The potential applied at the second
          barrier gate.
       plunger_gate_2_bias (float): The potential applied at the second
          plunger gate.
       barrier_gate_3_bias (float): The potential applied at the third
          barrier gate.

    Returns:
       Device: The Device object for the double quantum dot structure.

    """

    # Define the device object
    dvc = Device(mesh, conf_carriers="e")
    dvc.set_temperature(0.1)

    # Create the regions
    dvc.new_region("oxide", mt.SiO2)
    dvc.new_region("oxide_dot", mt.SiO2)
    dvc.new_region("gate_oxide", mt.HfO2)
    dvc.new_region("gate_oxide_dot", mt.HfO2)
    dvc.new_region("buried_oxide", mt.SiO2)
    dvc.new_region("buried_oxide_dot", mt.SiO2)
    dvc.new_region("channel", mt.Si)
    dvc.new_region("channel_dot", mt.Si)
    dvc.new_region("source", mt.Si, ndoping=1e20 * 1e6)
    dvc.new_region("drain", mt.Si, ndoping=1e20 * 1e6)

    # Set up boundary conditions
    Ew = mt.Si.Eg / 2 + mt.Si.chi  # Midgap
    dvc.new_gate_bnd("barrier_gate_1_bnd", barrier_gate_1_bias, Ew)
    dvc.new_gate_bnd("plunger_gate_1_bnd", plunger_gate_1_bias, Ew)
    dvc.new_gate_bnd("barrier_gate_2_bnd", barrier_gate_2_bias, Ew)
    dvc.new_gate_bnd("plunger_gate_2_bnd", plunger_gate_2_bias, Ew)
    dvc.new_gate_bnd("barrier_gate_3_bnd", barrier_gate_3_bias, Ew)
    dvc.new_ohmic_bnd("source_bnd")
    dvc.new_ohmic_bnd("drain_bnd")
    dvc.new_frozen_bnd(
        "back_gate_bnd", back_gate_bias, mt.Si, 1e15 * 1e6, "n", 46 * 1e-3 * ct.e
    )

    # Create the double quantum dot region
    dot_region_list = ["oxide_dot", "gate_oxide_dot", "buried_oxide_dot", "channel_dot"]
    dvc.set_dot_region(dot_region_list)

    return dvc
