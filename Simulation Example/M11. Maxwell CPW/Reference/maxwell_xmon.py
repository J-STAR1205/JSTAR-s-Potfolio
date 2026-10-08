__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

"""
Maxwell eigenmodes for an Xmon with the Josephson junction represented as a linear
inductor.
The layout and dimensions for this example were kindly provided by Christopher Xu from
Red Blue Quantum.
See the following article for more details on the operation of an Xmon qubit coupled to
a readout line, XY control line, and a quantum bus resonator.
    Barends, Rami, et al. "Coherent Josephson qubit suitable for scalable quantum
    integrated circuits." Phys. Rev. Lett., 111.8 (2013): 080502.
"""
from pathlib import Path
import os
from time import time
import numpy as np

# Import relevant modules of QTCAD.
from qtcad.device.maxwell_eigenmode import Solver
from qtcad.device.maxwell_eigenmode import SolverParams
from qtcad.device.device import Device
from qtcad.device.mesh3d import Mesh
from qtcad.device import materials as mt
from qtcad.device import constants as ct

# Scale in the Gmsh files. That is, Gmsh file coordinates are in μm.
scale = 1e-6

# Directories and file paths.
script_dir = Path(__file__).parent.resolve()
# For mesh and raw geometry files.
input_dir = script_dir / "meshes"
# For results.
result_dir = script_dir / "output" / Path(__file__).stem
# Mesh file.
fpath_mesh = input_dir / "xmon.msh"
# Raw geometry file (needed for adaptive meshing).
# NOTE: pointed at an ASCII-only path instead of input_dir because gmsh's
# XAO/XML reader fails to open files when the path contains non-ASCII
# (Korean) characters, as this project directory does. See tunnel_coupling_1.py
# for the same workaround. Content is identical to input_dir/xmon.xao.
fpath_xao = Path(
    r"C:\Users\norma\AppData\Local\Temp\claude\C--Users-norma-Desktop------QTCAD-Simulation\1c089a01-c171-407b-9be8-42b4f3db9e1e\scratchpad\ascii_geo\xmon.xao"
)

# Check if the mesh and raw geometry files exist.
if not os.path.isfile(fpath_mesh) or not os.path.isfile(fpath_xao):
    raise Exception(
        "Please run %s/xmon.py to generate the mesh and raw geometry files."
        % (input_dir)
    )

######################################################################################
# Setup the device.
######################################################################################
# Parse the mesh and initialize the device.
mesh = Mesh(scale, fpath_mesh)
dvc = Device(mesh)

material_sub = mt.Si
material_air = mt.vacuum
# Assign media to regions.
dvc.new_region("substrate", material_sub)
dvc.new_region("air", material_air)
# Assign perfect electric conductor boundary condition to all conductors.
dvc.new_pec_bnd("gnd")
dvc.new_pec_bnd("xmon_cross")
dvc.new_pec_bnd("xy_ctrl")
dvc.new_pec_bnd("readout")
dvc.new_pec_bnd("qbus")


# Calculate inductance.
#
# For the expression of the inductance as a function of the flux, see
#   Alexandre Blais, Arne L. Grimsmo, S. M. Girvin, and Andreas Wallraff.
#   Circuit quantum electrodynamics. Rev. Mod. Phys., 93:025005, May 2021.
EJ = 23e9 * ct.h
Phi0 = ct.h / (2 * ct.e)
inductance = Phi0**2 / ((2 * np.pi) ** 2 * EJ)

# Add the inductor.
dvc.new_inductor("jj", inductance, dir="y", length=24e-6, width=24e-6)

######################################################################################
# Setup the solver.
######################################################################################
params = SolverParams()
# Number of modes to find.
params.num_modes = 1
# Adaptive meshing (relative) tolerance on the frequency.
params.tol_rel = 0.03
# Directory where output files will be stored
params.output_dir = result_dir

######################################################################################
# Run the solver (results are stored in the folder specified by params.output_dir).
######################################################################################
slv = Solver(dvc, params, geo_file=fpath_xao)
t0 = time()
slv.solve()
dt = time() - t0
print("Solution completed in %.2f s" % dt)


######################################################################################
# Extract additional information from the final results.
######################################################################################

# From the cap_xmon.py tutorial.
capacitance = 1.014e-13

freq_lc = 1 / np.sqrt(inductance * capacitance) / (2 * np.pi)
freq_lc_ghz = freq_lc / 1e9
print(f"Frequency of the LC resonator (analytical): {freq_lc_ghz:.3f} GHz")

freq_ghz = dvc.maxwell_freqs[0] / 1e9
print(f"Frequency of the fundamental mode (computed): {freq_ghz:.3f} GHz")

# Blais et al. (2021).
EC = ct.e**2 / (2 * capacitance)
freq_qbit = (np.sqrt(8 * EC * EJ) - EC) / ct.h
freq_qbit_ghz = freq_qbit / 1e9
print(f"Frequency of the qubit (analytical): {freq_qbit_ghz:.3f} GHz")

from qtcad.device.epr import EPRAnalysis
from qtcad.device.epr import JunctionSpec, DielectricSpec, SurfaceDielectricSpec

epr_junctions = [JunctionSpec("jj")]

spec_dielectrics = [
    DielectricSpec("substrate", loss_tangent=5e-7),
]

spec_ms = SurfaceDielectricSpec(
    "xmon_cross",
    thickness=3e-9,
    relative_permittivity=10.0,
    loss_tangent=7e-4,
    interface_type="ms",
)
spec_ma = SurfaceDielectricSpec(
    "xmon_cross",
    thickness=3e-9,
    relative_permittivity=10.0,
    loss_tangent=4e-3,
    interface_type="ma",
)
spec_sa = SurfaceDielectricSpec(
    "substrate_top",
    thickness=3e-9,
    relative_permittivity=4,
    loss_tangent=6e-4,
    interface_type="sa",
)

######################################################################################
# Solve EPR analysis.
######################################################################################
# The EPR analysis reads the solved Maxwell fields and the current through the
# inductor boundary named "jj", then converts the linear EPR data into
# first-order nonlinear quantities such as the dressed frequency, zero-point
# fluctuations and anharmonicity.
epr = EPRAnalysis(
    dvc,
    epr_junctions,
    dielectrics=spec_dielectrics,
    surface_dielectrics=[spec_ms, spec_ma, spec_sa],
).solve()

######################################################################################
# EPR results and comparisons with simple analytical estimates.
######################################################################################
# Calculate relaxation time T1 from the total quality factor.
freq_dressed = epr.dressed_frequencies[0]
q_total = epr.mode_quality_factors[0]
t1 = q_total / (2 * np.pi * freq_dressed)

rows = [
    ("Frequency of the LC resonator (analytical)", f"{freq_lc / 1e9:.3f} GHz"),
    ("Frequency of the qubit (analytical)", f"{freq_qbit / 1e9:.3f} GHz"),
    ("Bare Maxwell frequency (computed)", f"{epr.bare_frequencies[0] / 1e9:.3f} GHz"),
    ("Dressed qubit frequency (EPR)", f"{freq_dressed / 1e9:.3f} GHz"),
    ("Junction participation", f"{epr.participation[0, 0]:.6f}"),
    ("Reduced flux ZPF", f"{epr.phi_zpf[0, 0]:.6f}"),
    ("Anharmonicity", f"{epr.anharmonicity[0] / 1e6:.3f} MHz"),
    ("Bulk substrate quality factor", f"{epr.dielectric_quality_factors[0, 0]:.3e}"),
    (
        "Surface MS quality factor",
        f"{epr.surface_dielectric_quality_factors[0, 0]:.3e}",
    ),
    (
        "Surface MA quality factor",
        f"{epr.surface_dielectric_quality_factors[0, 1]:.3e}",
    ),
    (
        "Surface SA quality factor",
        f"{epr.surface_dielectric_quality_factors[0, 2]:.3e}",
    ),
    ("Total quality factor", f"{q_total:.3e}"),
    ("Relaxation time T1", f"{t1 * 1e6:.3f} μs"),
]

# Display data as a table.
header = f"{'Quantity':<45} | {'Value':>15}"
print("\n" + header)
print("-" * len(header))
for label, value in rows:
    print(f"{label:<45} | {value:>15}")
