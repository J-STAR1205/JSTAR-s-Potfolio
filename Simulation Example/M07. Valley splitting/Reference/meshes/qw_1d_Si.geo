SetFactory("OpenCASCADE");

// Ignore previously user-defined (system) Gmsh options (files ~/.gmshrc and ~/.gmsh-options).
// This step ensures that the characteristic element lengths specified in the geometry file are used for meshing.
Delete Options;
General.Terminal = 1;
General.Verbosity = 5;

// Characterisitc length of the mesh
h           =   0.01;

// Define the points
Point(1)    =   {-13, 0, 0, h};
Point(2)    =   {0, 0, 0, h};   // Si/SiO2 interface
Point(3)    =   {3, 0, 0, h};
// Lines between points
Line(1)     =   {1, 2};
Line(2)     =   {2, 3};
// Physical lines - used to set material properties of regions
Physical Line("Domain", 1)      = {1};
Physical Line("Barrier", 2)     = {2};

Mesh 1;

Save "Si_well.msh";