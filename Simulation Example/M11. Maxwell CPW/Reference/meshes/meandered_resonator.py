__copyright__ = "2025, Nanoacademic Technologies Inc."

import gmsh
from pathlib import Path
import meandered_resonator_helper as helper
import math
import numpy as np

# Determine the directory the script is in.
script_dir = Path(__file__).parent.resolve()

########################################################################################
# Geometry parameters (length units: μm).
########################################################################################
# Length of the resonator.
length = 5950
# Thickness of the substrate.
sub_height = 500
# Thickness of the air layer above the substrate.
air_height = 500
# Position of the shorted end.
x0 = 500
y0 = 500
# Size of the chip in x and y dimensions.
lenx = 2500
leny = 1500
# Width of the central strip.
strip_width = 15
# Gap between the strip and the ground.
gap_width = 9
# Radius of the bends of the meandered resonator.
bend_rad = 90
# Number of straight line sections in the meandered resonator.
num_lines = 4

########################################################################################
# Description of each part of the resonator.
########################################################################################
# Number of bends.
n_bends = num_lines - 1
# Lenght of each 180-degree bend.
len_180_bends = math.pi * bend_rad
# Length of straight segments.
len_each = (length - n_bends * len_180_bends) / num_lines

# Meander description (See `helper.add_meander` for details).
meander = []
for i in range(num_lines):
    if i >= 1:
        if i % 2 == 1:
            meander.append("ccw")
            meander.append("ccw")
        else:
            meander.append("cw")
            meander.append("cw")
    meander.append(len_each)

########################################################################################
# Conductor crossovers along the y direction.
# They are modelled as air bridges with a square cross section.
########################################################################################
# Number of crossovers for each straight extension.
crossover_n = 2
# Side length of the conductor crossover’s cross section.
crossover_side = strip_width / 3
# Distance between the crossover bases to the edge of the ground plane close to the CPW.
crossover_gap = gap_width

########################################################################################
# Mesh parameters.
########################################################################################
# Characteristic lengths. `h_min` should be small enough to allow proper meshing of the
# conductor crossovers.
h_min = strip_width
h_max = 2000
# Calculate mesh element sizes based on curvature.
mesh_elements_curvature = 20

########################################################################################
# Mesh generation.
########################################################################################
# Initialize Gmsh (ignoring previously user-defined/system Gmsh options) and
# add the model.
gmsh.initialize(readConfigFiles=False)
gmsh.option.setNumber("General.Terminal", 1)
gmsh.option.setNumber("General.Verbosity", 5)
gmsh.model.add(Path(__file__).stem)

# Add substrate and air entities.
sub_tag = gmsh.model.occ.add_box(0, 0, -sub_height, lenx, leny, sub_height)
air_tag = gmsh.model.occ.add_box(0, 0, 0, lenx, leny, air_height)

# Add the entity for the ground plane.
gnd_plane_tag = gmsh.model.occ.add_rectangle(0, 0, 0, lenx, leny)

# Add the entity for the ground plane.
strip_dimtags, straight_sections = helper.add_meander(
    meander, strip_width, bend_rad, x0, y0
)

# Add the entity for the gap in the ground plane.
cut_dimtags, _ = helper.add_meander(
    [*meander, gap_width], strip_width + 2 * gap_width, bend_rad, x0, y0
)

gmsh.model.occ.synchronize()

# Origin of the conductor crossover start point.
crossover_origins = list()
crossovers = list()
pipe_dimtags = list()
crossover_dimtags = list()
# Get the locations to add conductor crossovers along the y direction.
for straight_section in straight_sections:
    c_x0, c_y0, c_z0, c_Δx, c_Δy = straight_section
    crossover_z0 = 0
    # Get points evenly distributed along the section, but away from its edges.
    crossover_x0s = np.linspace(c_x0, c_x0 + c_Δx, crossover_n + 1, endpoint=False)[1:]
    # Position the starting point below (xy plane) the strip.
    crossover_y0 = c_y0 - gap_width - (crossover_gap + crossover_side / 2)
    crossover_origins.extend(
        [(crossover_x0, crossover_y0, crossover_z0) for crossover_x0 in crossover_x0s]
    )

# Add crossovers.
for origin in crossover_origins:
    crossover_x0, crossover_y0, crossover_z0 = origin

    # Target end point of the conductor crossover.
    crossover_x1 = crossover_x0
    crossover_y1 = (
        crossover_y0
        + strip_width
        + (1 + 2) * gap_width
        + (crossover_gap + crossover_side)
    )
    # Radius of the arc.
    crossover_arc_radius = (crossover_y1 - crossover_y0) / 2

    gnd_plane_tag, pipe_tag, crossover_surf_tags = helper.add_crossover_arc_on_plane(
        origin=origin,
        arc_radius=crossover_arc_radius,
        side_length=crossover_side,
        tag_plane=gnd_plane_tag,
    )
    crossover_dimtag = [(2, t) for t in crossover_surf_tags]
    crossover_dimtags.extend(crossover_dimtag)
    pipe_dimtags += [(3, pipe_tag)]

# Fragment all entities involved.
to_fragment = [
    (3, sub_tag),
    (3, air_tag),
    *pipe_dimtags,
    (2, gnd_plane_tag),
    *crossover_dimtags,
    *strip_dimtags,
    *cut_dimtags,
]

_, fragment_map = gmsh.model.occ.fragment(to_fragment, [])
gmsh.model.occ.synchronize()

# Use fragment_map output of fragment function to determine the tags after
# fragmentation.
_, sub_tag = fragment_map.pop(0)[0]
_, air_tag = fragment_map.pop(0)[0]
pipe_dimtags = helper.pop_dimtags(fragment_map, len(pipe_dimtags))
# Note that the ground plane still includes the cutout.
gnd_plane_dimtags = fragment_map.pop(0)
crossover_dimtags = helper.pop_dimtags(fragment_map, len(crossover_dimtags))
strip_dimtags = helper.pop_dimtags(fragment_map, len(strip_dimtags))
cut_dimtags = helper.pop_dimtags(fragment_map, len(cut_dimtags))

# Exclude cutout entities from the ground plane.
gnd_plane_dimtags = list(set(gnd_plane_dimtags).difference(cut_dimtags))

# Obtain tags of the outer boundary (envelope).
envelope_tags = helper.get_gmsh_tags(
    gmsh.model.get_boundary(gmsh.model.get_entities(3), oriented=False)
)

# Remove crossover volumes.
gmsh.model.occ.remove(pipe_dimtags)
gmsh.model.occ.synchronize()

composite_gnd_plane_dimtags = [
    *helper.get_gmsh_tags(gnd_plane_dimtags),
    *helper.get_gmsh_tags(crossover_dimtags),
]

# Assign physical groups to 3D regions and 2D boundaries.
gmsh.model.add_physical_group(3, [sub_tag], name="substrate")
gmsh.model.add_physical_group(3, [air_tag], name="air")
gmsh.model.add_physical_group(2, composite_gnd_plane_dimtags, name="gnd")
gmsh.model.add_physical_group(2, helper.get_gmsh_tags(strip_dimtags), name="strip")
gmsh.model.add_physical_group(2, envelope_tags, name="envelope")

gmsh.model.occ.synchronize()

# Set characteristic length and mesh algorithm.
gmsh.option.setNumber("Mesh.MeshSizeMin", h_min)
gmsh.option.setNumber("Mesh.MeshSizeMax", h_max)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", mesh_elements_curvature)
# Use the HXT algorithm.
gmsh.option.setNumber("Mesh.Algorithm3D", 10)

# Generate and write the mesh.
gmsh.model.mesh.generate(3)
for ext in [".msh4", ".xao"]:
    path = script_dir / (Path(__file__).stem + ext)
    gmsh.write(str(path))

# Finalize Gmsh.
gmsh.finalize()
