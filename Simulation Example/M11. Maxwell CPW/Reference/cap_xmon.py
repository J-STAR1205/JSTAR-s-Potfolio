__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

# Capacitance matrix for an XMON with coupling capacitors.
#
# The layout and dimensions for this example were kindly provided by Christopher
# Xu, Red Blue Quantum.
#
# See the following article for more details on the operation of an Xmon qubit
# coupled to a readout line, XY control line, and a quantum bus resonator.
#
#    Barends, Rami, et al. "Coherent Josephson qubit suitable for scalable quantum
#    integrated circuits." Physical review letters 111.8 (2013): 080502.

from pathlib import Path
from time import time
import os

# import relevant modules of QTCAD
import qtcad.device.capacitance as cap
from qtcad.device.mesh3d import Mesh
from qtcad.device import Device
from qtcad.device import materials as mt

# directories and file paths
script_dir = Path(__file__).parent.resolve()
input_dir = script_dir / "meshes"  # for mesh and raw geometry files
result_dir = script_dir / "output" / Path(__file__).stem  # for results
fpath_mesh = input_dir / "xmon.msh"  # mesh file
# raw geometry file (needed for adaptive meshing)
# NOTE: pointed at an ASCII-only path instead of input_dir because gmsh's
# XAO/XML reader fails to open files when the path contains non-ASCII
# (Korean) characters, as this project directory does. See tunnel_coupling_1.py
# for the same workaround. Content is identical to input_dir/xmon.xao.
fpath_xao = Path(
    r"C:\Users\norma\AppData\Local\Temp\claude\C--Users-norma-Desktop------QTCAD-Simulation\1c089a01-c171-407b-9be8-42b4f3db9e1e\scratchpad\ascii_geo\xmon.xao"
)

# code to check that the mesh and raw geometry files exist
if not os.path.isfile(fpath_mesh) or not os.path.isfile(fpath_xao):
    raise Exception(
        "Please run %s/xmon.py to generate the mesh and raw geometry files."
        % (input_dir)
    )


# scale in the Gmsh files
scale = 1e-6  # i.e. Gmsh files are in um

# parse the mesh
mesh = Mesh(scale, fpath_mesh)

# Create device from mesh.
dvc = Device(mesh)

# Create regions
dvc.new_region("substrate", mt.Si)
dvc.new_region("air", mt.vacuum)

# Set potential to zero on the surfaces of each conductor (including ground)
dvc.new_dirichlet_bnd("gnd", 0.0)
dvc.new_dirichlet_bnd("xmon_cross", 0.0)
dvc.new_dirichlet_bnd("xy_ctrl", 0.0)
dvc.new_dirichlet_bnd("readout", 0.0)
dvc.new_dirichlet_bnd("qbus", 0.0)

######################################################################################
# setup and run the solver
######################################################################################

params = cap.SolverParams()

# directory where results will be stored
params.output_dir = result_dir
# tolerances
params.tol_rel = 0.05
params.tol_abs = 0
# name
params.name = Path(__file__).stem

# Set up signal conductors
signal_conductors = ["xmon_cross", "xy_ctrl", "readout", "qbus"]


# initialize the solver
cap_solver = cap.Solver(dvc, signal_conductors, params, geo_file=fpath_xao)

# execute the solver algorithm
t0 = time()
c = cap_solver.solve()
dt = time() - t0

######################################################################################
# Display the final results (in addition to the data in convergence files)
######################################################################################
cap_solver.print_cap(c, decimals=3)
print("Time taken: %.1f s" % dt)
