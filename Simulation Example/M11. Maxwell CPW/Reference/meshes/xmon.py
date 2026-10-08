__copyright__ = "2024, Nanoacademic Technologies Inc."

"""
Mesh generation for an XMON with coupling capacitors.

The layout and dimensions for this example were kindly provided by Christopher Xu,
Red Blue Quantum.

See the following article for more details on the operation of an Xmon qubit
coupled to a readout line, XY control line, and a quantum bus resonator.

    Barends, Rami, et al. "Coherent Josephson qubit suitable for scalable quantum
    integrated circuits." Physical review letters 111.8 (2013): 080502.
"""

import gmsh
from pathlib import Path
from pathlib import Path

import math
import xmon_helper

#############################################################################
# Define geometry parameters (length units: um)
#############################################################################

# metal thickness
#   Set to 0 to model metals as sheets
#   Set to a positive value to model the thickness of conductors
metal_t = 0

# height of the air meshed above the metal layer
air_height = 500
# substrate height
sub_height = 500
# width of the strips in the Xmon arms
strip = 24
# widht of the gaps around the Xmon's conductor
gap = 24
# length of the Xmon arms
arm_length = 176
# Size of the substrate in the xy directions
sub_xy = 4000
# The distance between the edges of the etching for the transmon and the XY control line
xy_dist = 12
# The length of the XY control line included in the simulation
xy_len = 280
# The width of the strip of the XY control line
xy_strip = 10
# The gap etched out around the XY control line
xy_gap = 6
# Two of the gap widths, depending on the side, etched out around the XY control line
xy_gap1 = 9  # the side closest to the Xmon
xy_gap2 = 11  # the side farthest from the Xmon
# The distance between the edges of the etching for the transmon and the readout line
r_dist = 4
# The lenghth of the readout resonator/quantum bus included in the simulation
r_len_line = 200
# The extension of the fork in the two large coupling capacitors
r_len_fork = 66
# Widths (depending on the location) of the strip of the two large coupling capacitors
r_strip_line = 15  # for the resonator lines
r_strip_fork = 23  # for the forks
r_strip_join = 20  # for the part of the capacitors that join the two parts of each fork
# The gaps etched out around the two large coupling capacitors
r_gap_outer = 9  # outer side (except the gap described by r_gap_outer_fork)
r_gap_inner = 4  # for the side adjacent to the Xmon
r_gap_outer_fork = 13  # for the outer side of the fork

# characteristic length
h = 500

#############################################################################
# directories and file paths
#############################################################################
directory = Path(__file__).parent.resolve()
msh_file = directory / "xmon.msh"  # path to mesh
geo_file = directory / "xmon.xao"  # path to raw geometry file

#############################################################################
# generate the mesh
#############################################################################
dim_metal = 3 if metal_t else 2

# Initialize Gmsh (ignoring previously user-defined/system Gmsh options) and
# add the model.
gmsh.initialize(readConfigFiles=False)
gmsh.option.setNumber("General.Terminal", 1)
gmsh.option.setNumber("General.Verbosity", 5)
gmsh.model.add(Path(__file__).stem)

# half of the Xmon size including the etching
xmon_halfspan = strip / 2 + arm_length + gap

# add substrate and air
sub = gmsh.model.occ.add_box(
    -sub_xy / 2, -sub_xy / 2, -sub_height, sub_xy, sub_xy, sub_height
)
air = gmsh.model.occ.add_box(-sub_xy / 2, -sub_xy / 2, 0, sub_xy, sub_xy, air_height)
# add the ground plane
if metal_t:
    # thick metal: add as a box of thickness metal_t
    gnd = gmsh.model.occ.add_box(-sub_xy / 2, -sub_xy / 2, 0, sub_xy, sub_xy, metal_t)
else:
    # metal as a sheet
    gnd = gmsh.model.occ.add_rectangle(-sub_xy / 2, -sub_xy / 2, 0, sub_xy, sub_xy)

# Xmon cross #################################################################
# Create a Gmsh entity for the metal of the Xmon cross
xmon_cross = xmon_helper.add_cross(strip, arm_length, metal_t)
# Create a Gmsh entity for the cutout from the ground plane for the Xmon cross
xmon_cutout = xmon_helper.add_cross(strip + 2 * gap, arm_length, metal_t)
# Cut out the xmon_cutout from the ground (xmon_cutout entity deleted by default)
gnd = gmsh.model.occ.cut([(dim_metal, gnd)], [(dim_metal, xmon_cutout)])[0][0][1]

gmsh.model.occ.synchronize()

# XY control #################################################################

# Create a Gmsh entithy for the XY control metal
xy = xmon_helper.add_strip(
    x=strip / 2 + arm_length + gap + xy_dist + xy_gap1,
    y=-xy_strip / 2,
    dx=xy_len,
    dy=xy_strip,
    metal_t=metal_t,
)
# Create a Gmsh entity for the XY control cutout
xy_cutout = xmon_helper.add_strip(
    x=strip / 2 + arm_length + gap + xy_dist,
    y=-xy_strip / 2 - xy_gap,
    dx=xy_len + xy_gap1 + xy_gap2,
    dy=xy_strip + 2 * xy_gap,
    metal_t=metal_t,
)
# Cut out the xy_cutout from the ground
gnd = gmsh.model.occ.cut([(dim_metal, gnd)], [(dim_metal, xy_cutout)])[0][0][1]
gmsh.model.occ.synchronize()

# Readout line and coupling capacitor ########################################

# half of the width of Xmon strip, including the gap
xmon_halfwidth = strip / 2 + gap

# lowest x value of the metal of the readout line
pos_x_minus = -xmon_halfspan - r_dist - r_gap_inner - r_strip_join - r_len_line
# Create a Gmsh entity for the metal of the readout line
r_line = xmon_helper.add_fork(
    pos_x_minus=pos_x_minus,
    line_dx=r_len_line,
    line_dy=r_strip_line,
    fork_outer_dx=r_strip_join + r_len_fork,
    fork_outer_dy=2 * (xmon_halfwidth + r_dist + r_gap_inner + r_strip_fork),
    fork_inner_dx=r_len_fork,
    fork_inner_dy=2 * (xmon_halfwidth + r_dist + r_gap_inner),
    metal_t=metal_t,
)
# Create a Gmsh entity for readout line cutout
r_cutout = xmon_helper.add_fork(
    pos_x_minus=pos_x_minus - r_gap_outer,
    line_dx=r_len_line,
    line_dy=r_strip_line + 2 * r_gap_outer,
    fork_outer_dx=r_strip_join + r_len_fork + 2 * r_gap_outer,
    fork_outer_dy=2
    * (xmon_halfwidth + r_dist + r_gap_inner + r_strip_fork + r_gap_outer_fork),
    fork_inner_dx=r_len_fork + r_gap_outer - r_gap_inner,
    fork_inner_dy=2 * (xmon_halfwidth + r_dist),
    metal_t=metal_t,
)
# Cut out the r_cutout from the ground
gnd = gmsh.model.occ.cut([(dim_metal, gnd)], [(dim_metal, r_cutout)])[0][0][1]

# Quantum bus and coupling capacitor #########################################

# Create a Gmsh entity for the quantum bus metal
qbus = xmon_helper.add_fork(
    pos_x_minus=pos_x_minus,
    line_dx=r_len_line,
    line_dy=r_strip_line,
    fork_outer_dx=r_strip_join + r_len_fork,
    fork_outer_dy=2 * (xmon_halfwidth + r_dist + r_gap_inner + r_strip_fork),
    fork_inner_dx=r_len_fork,
    fork_inner_dy=2 * (xmon_halfwidth + r_dist + r_gap_inner),
    metal_t=metal_t,
)

# Create a Gmsh entity for the quantum bus cutout
qbus_cutout = xmon_helper.add_fork(
    pos_x_minus=pos_x_minus - r_gap_outer,
    line_dx=r_len_line,
    line_dy=r_strip_line + 2 * r_gap_outer,
    fork_outer_dx=r_strip_join + r_len_fork + 2 * r_gap_outer,
    fork_outer_dy=2
    * (xmon_halfwidth + r_dist + r_gap_inner + r_strip_fork + r_gap_outer_fork),
    fork_inner_dx=r_len_fork + r_gap_outer - r_gap_inner,
    fork_inner_dy=2 * (xmon_halfwidth + r_dist),
    metal_t=metal_t,
)

# Orient the quantum bus and its cutout correctly
gmsh.model.occ.rotate(
    [(dim_metal, qbus), (dim_metal, qbus_cutout)], 0, 0, 0, 0, 0, 1, math.pi / 2
)

# Cut qbus_cutout from the ground plane
gnd = gmsh.model.occ.cut([(dim_metal, gnd)], [(dim_metal, qbus_cutout)])[0][0][1]

gmsh.model.occ.synchronize()

# If the metal has thickness, cut out the metal volume from air ##############

if metal_t:
    air = gmsh.model.occ.cut(
        [(3, air)],
        [(3, gnd), (3, xy), (3, r_line), (3, qbus), (3, xmon_cross)],
        removeTool=False,
    )[0][0][1]
    gmsh.model.occ.synchronize()

# Josephson junction #########################################################

jj = gmsh.model.occ.add_rectangle(-strip / 2, arm_length + strip / 2, 0, strip, gap)
gmsh.model.occ.synchronize()

# Important step: fragment all entities involved #############################

# Call fragment function
_, map = gmsh.model.occ.fragment(
    [
        (3, sub),  # substrate
        (3, air),  # air
        (dim_metal, gnd),  # ground
        (dim_metal, xmon_cross),  # Xmon cross metal
        (dim_metal, xy),  # XY control metal
        (dim_metal, r_line),  # readout line and its coupling capacitor
        (dim_metal, qbus),  # quantum bus and its coupling capacitor
        (2, jj),  # Josephson junction
    ],
    [],
)
gmsh.model.occ.synchronize()

# Use the map output of fragment function to determine the tags after fragmentation
_, sub = map.pop(0)[0]  # get new substrate tag
_, air = map.pop(0)[0]  # get new air tag
_, gnd = map.pop(0)[0]  # get new ground tag
_, xmon_cross = map.pop(0)[0]  # get new tag for Xmon cross metal
_, xy = map.pop(0)[0]  # get new tag for new XY control line metal
_, r_line = map.pop(0)[0]  # get new tag for new readout line metal
_, qbus = map.pop(0)[0]  # get new tag for quantum bus metal
_, jj = map.pop(0)[0]  # get new tag for Josephson junction

# Define physical groups #####################################################

# substrate and air physical groups
gmsh.model.add_physical_group(3, [sub], name="substrate")
gmsh.model.add_physical_group(3, [air], name="air")

# Add physical group for the 2D surface of the substrate surrounding the Xmon.
delta_z = 1e-4
all_z0_surfs = gmsh.model.get_entities_in_bounding_box(
    -sub_xy / 2 - delta_z,
    -sub_xy / 2 - delta_z,
    -delta_z,
    sub_xy / 2 + delta_z,
    sub_xy / 2 + delta_z,
    delta_z,
    dim=2,
)
all_z0_tags = [e[1] for e in all_z0_surfs]
metal_tags = {gnd, xmon_cross, xy, r_line, qbus, jj}
sa_tags = [t for t in all_z0_tags if t not in metal_tags]
gmsh.model.add_physical_group(2, sa_tags, name="substrate_top")

# Assign physical groups to *surfaces* of metals
xmon_helper.set_surf_phys_group((dim_metal, xmon_cross), name="xmon_cross")
xmon_helper.set_surf_phys_group((dim_metal, xy), name="xy_ctrl")
xmon_helper.set_surf_phys_group((dim_metal, r_line), name="readout")
xmon_helper.set_surf_phys_group((dim_metal, qbus), name="qbus")
xmon_helper.set_surf_phys_group((2, jj), name="jj")

# Special treatment of the ground phhysical group
if metal_t:
    # The code below ensures that the four outside faces of the ground plane are not
    # labeled as gnd.
    gnd_surf = gmsh.model.get_boundary([(3, gnd)], oriented=False)
    outer_box = gmsh.model.occ.add_box(
        -sub_xy / 2, -sub_xy / 2, -sub_height, sub_xy, sub_xy, sub_height + air_height
    )
    gmsh.model.occ.synchronize()
    outer_box_bnd = gmsh.model.get_boundary([(3, outer_box)], oriented=False)
    # cut out the part of the ground metal surface where it touches the outer box
    # of the domain
    gnd_surf = gmsh.model.occ.cut(gnd_surf, outer_box_bnd)[0]
    gnd_surf = [dimtag[1] for dimtag in gnd_surf]
    gmsh.model.occ.remove([(3, outer_box)])
    gmsh.model.occ.synchronize()
    gmsh.model.add_physical_group(2, gnd_surf, name="gnd")
else:
    gmsh.model.add_physical_group(2, [gnd], name="gnd")

# Set characteristic lenght, verbosity, and Gmsh algorithm ###################
gmsh.option.setNumber("Mesh.MeshSizeMin", h)
gmsh.option.setNumber("Mesh.MeshSizeMax", h)
gmsh.option.setNumber("General.Verbosity", 1)
gmsh.option.setNumber("Mesh.Algorithm3D", 10)

# Create and save the mesh ###################################################
gmsh.model.mesh.generate(3)

for path in [msh_file, geo_file]:
    gmsh.write(str(path))

gmsh.finalize()
