// Ignore previously user-defined (system) Gmsh options (files ~/.gmshrc and ~/.gmsh-options).
// This step ensures that the characteristic element lengths specified in the geometry file are used for meshing.
Delete Options;
General.Terminal = 1;
General.Verbosity = 5;

// Geometry parameters
cap_thick = 10;
barrier_thick = 20;
substrate_thick = 10;

// Mesh layers
cap_layers = 500;
barrier_layers = 500;
substrate_layers = 500;

// Extrusions
Point(1) = {0, 0, 0, 1.0};
Extrude {cap_thick, 0, 0} {
  Point{1}; Layers {cap_layers};
}

Extrude {barrier_thick, 0, 0} {
  Point{2}; Layers {barrier_layers};
}

Extrude {substrate_thick, 0, 0} {
  Point{3}; Layers {substrate_layers};
}

// Regions
Physical Curve("left_barrier", 7) = {1};
Physical Curve("well", 8) = {2};
Physical Curve("right_barrier", 9) = {3};

// Generate and save mesh
Mesh 1;
Save "band_alignment.msh";
