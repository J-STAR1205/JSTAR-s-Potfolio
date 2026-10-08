// FDSOI periodic structure
// default mesh size
// Remarks: length units are nanometers

// Use built-in meshing kernel
SetFactory("OpenCASCADE");

// Ignore previously user-defined (system) Gmsh options (files ~/.gmshrc and ~/.gmsh-options).
// This step ensures that the characteristic element lengths specified in the geometry file are used for meshing.
Delete Options;
General.Terminal = 1;
General.Verbosity = 5;

// x axis dimensions
domain_width = 60;
channel_width = 50;

Mesh.CharacteristicLengthMin = 2;
Mesh.CharacteristicLengthMax = 5;

// y axis dimensions
gap_len_1 = 5;    // Length of the gap between edge and barrier 1
gap_len_2 = 5;    // Length of the gap between barrier 1 and plunger 1
gap_len_3 = 5;    // Length of the gap between plunger 1 and barrier2
gap_len_4 = .1;    // Length of the gap between barrier 2 and edge
clavier_gate_len = 30;    // Clavier gate length

// z axis dimensions
gate_oxide_thick = 10;
film_thick = 30;
box_thick = 30;

// Rectangle for the simulation domain
channel_len = gap_len_1 + gap_len_2 + gap_len_3 + gap_len_4 + 3*clavier_gate_len;
Rectangle(1) = {-domain_width/2, -channel_len/2, 0, domain_width, channel_len, 0}; //simulation domain

// Rectangle for the channel
Rectangle(2) = {-channel_width/2, -channel_len/2, 0, channel_width, channel_len, 0}; // channel area

// Rectangles for the top gates
gap_list = {gap_len_1, gap_len_2, gap_len_3, gap_len_4, 0};
temp = -channel_len/2+ gap_len_1;
For i In {3:5}
    Rectangle(i) = {-domain_width/2, temp  , 0, domain_width, clavier_gate_len, 0};
    temp = temp + gap_list[i-3] + clavier_gate_len; 
EndFor

// Boolean fragments
BooleanFragments{ Surface{1:5}; Delete; }{ }

// Extrude upwards to form the gate oxide
// Extrude function has an output that gives all the bottom surfaces
Extrude {0, 0, gate_oxide_thick} {Surface{1:21};}

// Physical surfaces for the gate boundaries
Physical Surface("A0") = {34, 41, 48};
Physical Surface("A1") = {55, 61, 68};
Physical Surface("A2") = {75, 81, 88};

// Physical volumes for the gate oxide
Physical Volume("gate_oxide") = {1:21};

// Extrude downwards to form the film
Extrude {0, 0, -film_thick} {Surface{1:21};}

// Physical volumes for the channel and oxide
Physical Volume("channel") = {23,26,29,32,35,38,41};
Physical Volume("oxide") = {22,24,27,30,33,36,39,25,28,31,34,37,40,42};

// Extrude downwards to form the buried oxide
Extrude {0, 0, -box_thick} {
    Surface{99}; Surface{103}; Surface{107}; Surface{111};
    Surface{114}; Surface{118}; Surface{121}; Surface{124}; 
    Surface{128}; Surface{131}; Surface{134}; Surface{138}; 
    Surface{141}; Surface{144}; Surface{148}; Surface{151}; 
    Surface{154}; Surface{158}; Surface{161}; Surface{164}; 
    Surface{167};}

// Physical volume for the buried oxide
Physical Volume("buried_oxide") = {43:63};

// Physical surface for the back gate
Physical Surface("back_gate_bnd") = {172,176,180,184,187,191,194,197,201
    ,204,207,211,214,217,221,224,227,231,234,237,240};

// Physical surfaces for the Periodic boundaries
// Define surface ID lists
right_bnd[] = {83,89,92,156,162,165,229,235,238}; // Right boundary (follower)
left_bnd[]  = {25,29,37,98,102,110,171,175,183};  // Left boundary (main)

// Apply periodicity for mesh
Periodic Surface{right_bnd[]} = {left_bnd[]} Translate {0, channel_len, 0};

// Define Physical Surfaces for the periodic boundaries
Physical Surface("right_bnd") = {right_bnd[]};
Physical Surface("left_bnd")  = {left_bnd[]};

// Mesh
Mesh 1;
Mesh 2;
Mesh 3;

//Save
Save "periodic_fdsoi.msh";
Save "periodic_fdsoi.xao";