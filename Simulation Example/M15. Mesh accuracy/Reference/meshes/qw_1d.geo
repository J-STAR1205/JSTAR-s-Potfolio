// Ignore previously user-defined (system) Gmsh options (files ~/.gmshrc and ~/.gmsh-options).
// This step ensures that the characteristic element lengths specified in the geometry file are used for meshing.
Delete Options;
General.Terminal = 1;
General.Verbosity = 5;

// Geometry parameters
well_width = 2;

// Mesh parameters
num_elems = 250;

// Create structure
Point(1) = {0, 0, 0, 1.0};
Extrude {well_width, 0, 0} {
  Point{1}; Layers {num_elems}; 
}
Physical Point("left_barrier", 1) = {1};
Physical Point("right_barrier", 2) = {2};
Physical Curve("well", 6) = {1};

Mesh 1;
Save "quantum_well.msh";