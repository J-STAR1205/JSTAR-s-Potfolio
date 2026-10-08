import gmsh
import numpy as np


def add_crossover_arc_on_plane(
    origin: tuple[float, float, float],
    arc_radius: float,
    side_length: float,
    tag_plane: int,
    n_discret_points: int = 100,
) -> tuple[int, int, list[int, ...]]:
    x0, y0, z0 = origin
    # Storage for the points to make up the wire.
    points4spline = list()
    dimtags_to_remove_recursive = list()
    # Add points to create a spline and then a wire.
    for i in range(n_discret_points):
        # (n_discret_points - 1) to get the pipe level to the plane.
        theta = i * np.pi / (n_discret_points - 1)
        # Coordinates.
        _x = x0
        _y = (y0 + arc_radius) - arc_radius * np.cos(theta)
        _z = z0 + arc_radius * np.sin(theta)
        # ‘Punch’ the last point to avoid fragmentation issues.
        if i == n_discret_points - 1:
            _z = z0 - 1e-4

        tag_point = gmsh.model.occ.addPoint(_x, _y, _z)
        points4spline.append(tag_point)
        dimtags_to_remove_recursive.append((0, tag_point))
    tag_spline = gmsh.model.occ.addSpline(points4spline)
    dimtags_to_remove_recursive.append((1, tag_spline))
    tag_wire = gmsh.model.occ.addWire([tag_spline])

    # Insert the rectangle to form the pipe along the wire.
    tag_disk1 = gmsh.model.occ.addRectangle(x0, y0, z0, side_length, side_length)
    tag_pipe = gmsh.model.occ.addPipe([(2, tag_disk1)], tag_wire)[0][-1]
    gmsh.model.occ.synchronize()

    # Fragment everything.
    _, fragment_map = gmsh.model.occ.fragment(
        [(3, tag_pipe), (2, tag_disk1), (2, tag_plane)], []
    )
    pipe_fragments = list(
        sorted(fragment_map[0], key=lambda dimtag: -gmsh.model.occ.getMass(*dimtag))
    )
    # Remove everything but the largestbit.
    gmsh.model.occ.remove(pipe_fragments[1:], recursive=True)
    gmsh.model.occ.synchronize()
    # Recover tags.
    _, tag_pipe = pipe_fragments[0]
    fragment_map.pop(0)
    _, tag_disk1 = fragment_map.pop(0)[0]
    _, tag_plane = fragment_map.pop(0)[0]
    gmsh.model.occ.synchronize()

    # Get the pipe’s surface.
    _, surf_list_pipe = gmsh.model.getAdjacencies(3, tag_pipe)

    # Remove unwanted objects. The pipe’s volume should remain as we need to
    # make everything conformal afterwards, including any surrounding
    # dielectric media.
    gmsh.model.occ.remove(dimtags_to_remove_recursive, recursive=True)

    return tag_plane, tag_pipe, surf_list_pipe


def pop_dimtags(
    fragment_map: list[list[tuple[int, int]]], num_dimtags: int
) -> list[tuple[int, int]]:
    """Interprets the map that is the second output in gmsh.model.occ.fragment to find
    the correspondence between a collection of old entities and the new entities.
    Removes ("pops") the corresponding information from the map.

    Args:
        fragment_map: list where each element describes the result gmsh.model.occ.fragment
            had on each entity passed to that function. Each element of fragment_map is
            itself a list of tuples of new dimtags corresponding to the original entity.
        num_dimtags: how many entries to extract from fragment_map.

    Returns:
        list of new (dim, tag) tuples that correspond to the first num_dimtags elements of
        fragment_map.
    """
    dimtags = []
    for _ in range(num_dimtags):
        dimtags.extend(fragment_map.pop(0))
    return dimtags


def cw(dir: np.ndarray[float]) -> np.ndarray[float]:
    """Find direction clockwise from the given direction.

    Args:
        dir: 2-element 1D array representing a vector in the xy plane.
            The vector specifies the direction in which the resonator
            is being currently extended.

    Returns:
        2-element 1D array representing a vector pointing in the
        direction in which the resonator would continue if it were to
        bend clockwise by 90 degrees.
    """
    return np.cross([*dir, 0], [0, 0, 1])[:2]


def ccw(dir):
    """Find direction counterclockwise from the given direction.

    Args:
        dir: 2-element 1D array representing a vector in the xy plane.
            The vector specifies the direction in which the resonator
            is being currently extended.

    Returns:
        2-element 1D array representing a vector pointing in the
        direction in which the resonator would continue if it were to
        bend counterclockwise by 90 degrees.
    """
    return -cw(dir)


def add_bend(
    xc: float,
    yc: float,
    quadrant: np.ndarray[float],
    width: float,
    r: float,
):
    """Add a bend to the resonator.

    Args:
        xc, yc: current position before adding the bend.
        quadrant: a 2-element vector pointing at 45 degrees into the quadrant
            where the bend lies, when looking from the center point of the arc.
        width: the width of the line of the bend.
        r: the radius of the bend.
    """
    ro = r + width / 2
    ri = r - width / 2
    disk_outer = gmsh.model.occ.add_disk(xc=xc, yc=yc, zc=0, rx=ro, ry=ro)
    disk_inner = gmsh.model.occ.add_disk(xc=xc, yc=yc, zc=0, rx=ri, ry=ri)
    x1 = xc
    y1 = yc
    x2 = xc + 1.1 * ro * quadrant[0]
    y2 = yc + 1.1 * ro * quadrant[1]

    x = min(x1, x2)
    y = min(y1, y2)
    dx = max(x1, x2) - x
    dy = max(y1, y2) - y

    # arbitrary value of similar order of magnitude as others in the problem
    delta = max(abs(dx), abs(dy))

    intersect = gmsh.model.occ.add_box(x, y, -delta, dx, dy, 2 * delta)

    gmsh.model.occ.synchronize()

    _, map = gmsh.model.occ.intersect(
        [(2, disk_outer), (2, disk_inner)], [(3, intersect)]
    )

    gmsh.model.occ.synchronize()

    _, disk_inner = map[1][0]
    disk_outer = np.setdiff1d(get_gmsh_tags(map[0]), disk_inner)[0]

    disk = gmsh.model.occ.cut([(2, disk_outer)], [(2, disk_inner)])[0][0][1]
    gmsh.model.occ.synchronize()
    return disk


def add_meander(
    meander: list[str | float],
    width: float,
    bend_rad: float,
    x0: float,
    y0: float,
) -> tuple[list[tuple[int, int]], list[tuple[float, float, float, float, float]]]:
    """Adds a meandering line.

    Args:
        meander: a list of elements of the following format:
            -   string "cw" or "ccw" indicating 90-degree bends and their
                directions (clockwise or counterclockwise).
            -   scalar indicating a straight extension of the given length.
        width: the width of the line of the meander.
        bend_rad: the radius of the bends of the meander.
        x0, y0: the initial position of the meander in the xy plane.
    """
    cur_pos = np.array([x0, y0], dtype=float)
    cur_dir = np.array([1, 0], dtype=float)  # +x
    dir_fncs = dict(cw=cw, ccw=ccw)
    dimtags = []
    straight_sections = list()
    for val in meander:
        if isinstance(val, str):
            new_dir = dir_fncs[val](cur_dir)
            pos_c = cur_pos + new_dir * bend_rad
            tag = add_bend(pos_c[0], pos_c[1], cur_dir - new_dir, width, bend_rad)

            cur_pos += cur_dir * bend_rad + new_dir * bend_rad
            cur_dir = dir_fncs[val](cur_dir)
        elif isinstance(val, float) | isinstance(val, int):
            dx, dy = cur_dir * val + abs(cw(cur_dir)) * width
            x, y = cur_pos - abs(cw(cur_dir)) * width / 2
            coordinates = (x, y, 0)
            tag = gmsh.model.occ.add_rectangle(*coordinates, dx, dy)
            gmsh.model.occ.synchronize()
            straight_sections += [(*coordinates, dx, dy)]
            cur_pos += cur_dir * val
        else:
            raise Exception("Invalid type")
        dimtags.append((2, tag))

    gmsh.model.occ.synchronize()

    return dimtags, straight_sections


def get_gmsh_tags(dim_tags: list[tuple[int, int]]) -> list[int]:
    """From a list of (dim, tag) tuples, extract a list of tags.

    Args:
        dim_tags: the list of (dim, tag) tuples.

    Returns:
        list containing just the tag integers from the original list.
    """
    tags = []
    for dim_tag in dim_tags:
        tags.append(dim_tag[1])
    return tags
