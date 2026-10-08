__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import pathlib
import numpy as np
from qtcad.device import constants as ct
from qtcad.device import materials as mt
from qtcad.device.mesh3d import Mesh, SubMesh
from qtcad.device import Device
from qtcad.device import io

script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = str(script_dir / "meshes" / "multiscale.msh")
scaling = 1e-9
mesh = Mesh(scaling, path_mesh)

d = Device(mesh, conf_carriers="e")

qd = ["cap_qd", "well_qd", "buffer_qd"]
submesh = SubMesh(mesh, qd)

path_der_hdf5 = str(script_dir / "output" / "multiscale_der.hdf5")
der = io.load(path_der_hdf5, var_name="der", asarray=True)
print("raw der min/max/mean:", der.min(), der.max(), der.mean())
der_meV = der / ct.e * 1000.0  # meV/V

coords = submesh.loc_node_coords()
print("loc_node_coords type:", type(coords), "len:", len(coords))
for i, c in enumerate(coords):
    print(f"  item {i}: type={type(c)}, shape={getattr(c, 'shape', None)}")
print("der shape:", der.shape)

x_loc, y_loc, z_loc = coords  # (n_elem, 4) each, assumed
x = x_loc.reshape(-1)
y = y_loc.reshape(-1)
z = z_loc.reshape(-1)
vals = der_meV.reshape(-1)

r = np.sqrt(x**2 + y**2)

# slice just below SiGe spacer / Si QW interface, same as tutorial (z = -36 nm)
z0 = -36.0 * scaling
dz = 1.0 * scaling  # +/- 1 nm tolerance
in_slice = np.abs(z - z0) < dz

r_slice = r[in_slice]
v_slice = vals[in_slice]

# center: r < 5 nm
center_mask = r_slice < 5.0 * scaling
# periphery: r between 40 and 50 nm (near edge of QD region, within plunger gate radius ~50nm)
periph_mask = (r_slice > 40.0 * scaling) & (r_slice < 50.0 * scaling)

print(f"Slice points: {in_slice.sum()}")
print(f"Center (r<5nm) points: {center_mask.sum()}, mean derivative = {v_slice[center_mask].mean():.4f} meV/V")
print(f"Periphery (40<r<50nm) points: {periph_mask.sum()}, mean derivative = {v_slice[periph_mask].mean():.4f} meV/V")
print(f"Overall slice min/max: {v_slice.min():.4f} / {v_slice.max():.4f} meV/V")
