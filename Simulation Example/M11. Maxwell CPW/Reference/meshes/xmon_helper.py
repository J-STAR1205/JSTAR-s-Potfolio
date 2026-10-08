import gmsh


def add_strip(x: float, y: float, dx: float, dy: float, metal_t: float) -> int:
    """Adds a rectangle with the bottom at z = 0 height.

    Args:
        x, y: coordinates of the corner of the rectangle with lower values
            of x and y.
        dx, dy: x and y dimensions of the rectangle.
        metal_t: thickness of the metal.

    Returns:
        Tag of the rectangle.
    """
    if metal_t:
        tag = gmsh.model.occ.add_box(x, y, 0, dx, dy, metal_t)
    else:
        tag = gmsh.model.occ.add_rectangle(x, y, 0, dx, dy)
    return tag


def add_cross(width: float, armlen: float, metal_t: float) -> int:
    """Adds a cross with center at (0,0,0).

    Args:
        width: the width of the arms of the cross.
        armlen: the lenght of the arms of the cross.
        metal_t: thickness of the metal (can be 0).

    Returns:
        The tag of the cross.
    """
    length = width + 2 * armlen

    hor = add_strip(
        x=-length / 2,
        y=-width / 2,
        dx=length,
        dy=width,
        metal_t=metal_t,
    )
    ver = add_strip(
        x=-width / 2,
        y=-length / 2,
        dx=width,
        dy=length,
        metal_t=metal_t,
    )
    dim_cross = 3 if metal_t else 2

    gmsh.model.occ.synchronize()
    cross, _ = gmsh.model.occ.fuse([(dim_cross, hor)], [(dim_cross, ver)])
    gmsh.model.occ.synchronize()
    assert len(cross) == 1
    return cross[0][1]


def add_fork(
    pos_x_minus: float,
    line_dx: float,
    line_dy: float,
    fork_outer_dx: float,
    fork_outer_dy: float,
    fork_inner_dx: float,
    fork_inner_dy: float,
    metal_t: float,
) -> int:
    """Add a fork with given parameters.

    Args:
        pos_x_minus: the x position of the beginning of the line of the fork.
        line_dx, line_dy: x and y size of the line (line_dx being length).
        fork_outer_dx, fork_outer_dx: the size of the outer rectangle of the fork.
        fork_inner_dy, fork_inner_dx: the size of the cutout of the fork.
        metal_t: thickness of the metal (can be 0).

    Returns:
        The tag of the fork.
    """
    dim_metal = 3 if metal_t else 2

    # add the line
    line = add_strip(pos_x_minus, -line_dy / 2, line_dx, line_dy, metal_t)
    # add the outer fork
    fork_outer = add_strip(
        pos_x_minus + line_dx, -fork_outer_dy / 2, fork_outer_dx, fork_outer_dy, metal_t
    )
    # add the cutout of the fork
    fork_inner = add_strip(
        pos_x_minus + line_dx + fork_outer_dx - fork_inner_dx,
        -fork_inner_dy / 2,
        fork_inner_dx,
        fork_inner_dy,
        metal_t,
    )
    gmsh.model.occ.synchronize()
    # create the outer-inner shape of the fork
    fork, _ = gmsh.model.occ.cut([(dim_metal, fork_outer)], [(dim_metal, fork_inner)])
    gmsh.model.occ.synchronize()
    assert len(fork) == 1
    fork = fork[0][1]  # tag of the fork
    # create the entire fork with the line
    fork_line, _ = gmsh.model.occ.fuse([(dim_metal, fork)], [(dim_metal, line)])
    gmsh.model.occ.synchronize()
    assert len(fork_line) == 1
    fork_line = fork_line[0][1]  # tag of line + fork
    return fork_line


def set_surf_phys_group(dimtag: tuple[int, int], name: str):
    """Sets a physical group for a surface. Can take as input a 2D or
    3D entity. For a 3D entity, sets the physical group for the entities on
    its boundary.

    Args:
        dimtag: dimension and the tag of the 2D or 3D entity.
        name: the name to give to the physical group.
    """
    dim, tag = dimtag
    if dim == 2:
        gmsh.model.add_physical_group(2, [tag], name=name)
    else:
        assert dim == 3
        bnd = gmsh.model.get_boundary([dimtag], oriented=False)
        gmsh.model.add_physical_group(
            2, [dimtag_bnd[1] for dimtag_bnd in bnd], name=name
        )
