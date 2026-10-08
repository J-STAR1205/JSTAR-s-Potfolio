// A solid ball of diameter 1 centered at the origin at the center of a cube
// of side length 2

SetFactory("OpenCASCADE");

// Ignore previously user-defined (system) Gmsh options (files ~/.gmshrc and ~/.gmsh-options).
// This step ensures that the characteristic element lengths specified in the geometry file are used for meshing.
Delete Options;
General.Terminal = 1;
General.Verbosity = 5;

// Parameters
r = 0.5; // ball radius
L = 2.0; // box side length
lc = 0.1; // characteristic mesh length

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

// Create a sphere (solid ball)
Sphere(1) = {0, 0, 0, r};

// Create a cube around it
Box(2) = {-L/2, -L/2, -L/2, L, L, L};

// Boolean difference
ball[] = Duplicata{ Volume{1}; };
cube_shell[] = BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; };

// Physical volumes
Physical Volume("ball")       = { ball[] };
Physical Volume("cube_shell") = { cube_shell[] };

Mesh 1;
Mesh 2;
Mesh 3;
Save "spherical_dot.msh";
