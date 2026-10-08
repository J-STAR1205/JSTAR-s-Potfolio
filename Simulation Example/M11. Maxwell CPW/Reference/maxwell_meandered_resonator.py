__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

"""
Find the first two modes of a meandered quarter-wavelength coplanar waveguide resonator
inside a cavity.
"""
from pathlib import Path
import os
from time import time
import numpy as np

# Import the relevant modules of QTCAD.
from qtcad.device.maxwell_eigenmode import Solver
from qtcad.device.maxwell_eigenmode import SolverParams
from qtcad.device.device import Device
from qtcad.device.mesh3d import Mesh
from qtcad.device import materials as mt
from qtcad.device import constants as ct

# Scale in the Gmsh files.
scale = 1e-6

# The total length of the resonator in meters.
length = 5950 * scale

# Directories and file paths.
script_dir = Path(__file__).parent.resolve()
# For mesh and raw geometry files.
input_dir = script_dir / "meshes"
# For results.
result_dir = script_dir / "output" / Path(__file__).stem
# Mesh file.
fpath_mesh = input_dir / "meandered_resonator.msh4"
# Raw geometry file.
# NOTE: pointed at an ASCII-only path instead of input_dir because gmsh's
# XAO/XML reader fails to open files when the path contains non-ASCII
# (Korean) characters, as this project directory does. See tunnel_coupling_1.py
# for the same workaround. Content is identical to input_dir/meandered_resonator.xao.
fpath_xao = Path(
    r"C:\Users\norma\AppData\Local\Temp\claude\C--Users-norma-Desktop------QTCAD-Simulation\1c089a01-c171-407b-9be8-42b4f3db9e1e\scratchpad\ascii_geo\meandered_resonator.xao"
)

# Check if the mesh and raw geometry files exist.
if not os.path.isfile(fpath_mesh) or not os.path.isfile(fpath_xao):
    raise Exception(
        "Please run %s/meandered_resonator.py to generate the mesh and raw geometry files."
        % (input_dir)
    )


######################################################################################
# Set up the device.
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
dvc.new_pec_bnd("strip")
dvc.new_pec_bnd("envelope")

######################################################################################
# Set up the solver.
######################################################################################
params = SolverParams()
# Number of modes to find.
params.num_modes = 2
# Adaptive meshing (relative) tolerance on the frequency.
params.tol_rel = 0.05
# Directory where output files will be stored.
params.output_dir = result_dir

# Uncomment the line below in order to save all intermediate results at every iteration.
# params.save_intermediate_results = True

# Number of consecutive iterations that must agree within the tolerance thresholds.
# Uncomment the line below to change the default behaviour.
# Then, for each mode, the last result will be compared to 7 previous data points (i.e.
# 8 data points in total are compared for each mode).
# params.min_converged_iters = 7

# Do not terminate until the refined mesh contains at least min_nodes nodes.
params.min_nodes = 20000


######################################################################################
# Run the solver (results are stored in the folder specified by params.output_dir).
######################################################################################
slv = Solver(dvc, params, geo_file=fpath_xao)
t0 = time()
slv.solve()
dt = time() - t0
print("Solution completed in %.2f s" % dt)

frequencies_ghz = [f"{f / 1e9:.3f}" for f in dvc.maxwell_freqs]
print("Frequencies found (GHz): " + ", ".join(frequencies_ghz))

######################################################################################
# Extract additional information from the final results.
######################################################################################

# Find the theoretical value
#
# NOTE:
#   The effective permittivity is computed for a coplanar waveguide placed on an
#   infinitely thick and wide substrate.
#
#       Simons, Rainee N. Coplanar waveguide circuits, components, and systems.
#       John Wiley & Sons, 2001.
#
#   Effects of the finite substrate size or the cavity are neglected due to the large
#   distance between the resonator and the cavity (and, consequently, also the edges of
#   the substrate).
#   Please refer to Simons (2001) for details on the treatment of such effects.
#   It should also be noted that the theoretical value does not take into account the
#   meandered layout of the resonator.
eps_eff = (material_air.eps + material_sub.eps) / 2
wavelength = length * 4
vphase = 1 / np.sqrt(eps_eff * ct.mu0)
freq_theor = vphase / wavelength
print(f"Frequency of the fundamental mode (analytical): {freq_theor / 1e9:.3f} GHz")

# Display the frequency of the fundamental mode and compare to theoretical value.
err = abs(dvc.maxwell_freqs[0] - freq_theor) / freq_theor
print(
    "\nFrequency of the fundamental mode (computed): %.3f GHz (rel. err: %.3f)"
    % (dvc.maxwell_freqs[0] / 1e9, err)
)

# For a quarter wavelength resonator, the resonances should occur at odd multiples of
# its fundamental frequency. For the first two modes, their ratio should be equal to 3.
freq_ratio = dvc.maxwell_freqs[1] / dvc.maxwell_freqs[0]
freq_ratio_theor = 3
err_freq_ratio = abs(freq_ratio - freq_ratio_theor) / freq_ratio_theor
print(
    f"\nRatio between the first excited and fundamental mode (computed): {freq_ratio:.3f}"
    f"\nRelative error from the expected ratio of {freq_ratio_theor}: {err_freq_ratio:.3f}"
)

# Get the electric and magnetic fields.
e_field = dvc.get_e_field_maxwell()
b_field = dvc.get_b_field_maxwell()

# Fields in the fundamental mode.
e_field_0 = e_field[..., 0]
b_field_0 = b_field[..., 0]

# Find the ratio between the energy stored in the substrate and the energy in the air.
energy_sub = dvc.energy_e(e_field_0, "substrate") + dvc.energy_b(b_field_0, "substrate")
energy_air = dvc.energy_e(e_field_0, "air") + dvc.energy_b(b_field_0, "air")

print(
    "\nIn the fundamental mode, the ratio of energy in the substrate and in the air is"
    f" {energy_sub / energy_air:.1f}:1"
)
