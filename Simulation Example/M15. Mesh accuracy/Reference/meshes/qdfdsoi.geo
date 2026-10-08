// half a Double Quantum dot in Fully Depleted Silicon On Insulator  (QDFDSOI)
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
channel_width = 40;

// y axis dimensions
gap_len_1 = 5;    // Length of the gap between source and barrier 1
gap_len_2 = 5;    // Length of the gap between barrier 1 and dot 1
gap_len_3 = 5;    // Length of the gap between dot 1 and barrier 2
plunger_gate_len = 15;    // Plunger gate length
barrier_gate_len = 10;     // Barrier gate length
source_drain_len = 20;     // Source/drain length

// z axis dimensions
gate_oxide_thick = 2;
film_thick = 10;
box_thick = 10;

// Rectangle for the simulation domain
channel_len = gap_len_1 + gap_len_2 + gap_len_3+ 1.5*barrier_gate_len + plunger_gate_len;
domain_len = source_drain_len + channel_len;
Rectangle(1) = {-domain_width/2, -domain_len/2, 0, domain_width, domain_len, 0}; //simulation domain
Printf("%f",domain_len);
// Rectangles for the source, channel and drain
Rectangle(2) = {-channel_width/2, -domain_len/2, 0, channel_width, source_drain_len, 0}; // source area
Rectangle(3) = {-channel_width/2, -domain_len/2 +  source_drain_len, 0, channel_width, channel_len, 0}; // channel area

// Rectangles for the plunger gates and barrier gates
Rectangle(4) = {-domain_width/2, -domain_len/2 +  source_drain_len + gap_len_1  , 0, domain_width, barrier_gate_len, 0};
Rectangle(5) = {-domain_width/2, -domain_len/2 +  source_drain_len + gap_len_1 + barrier_gate_len + gap_len_2  , 0, domain_width, plunger_gate_len, 0};
Rectangle(6) = {-domain_width/2, -domain_len/2 +  source_drain_len + gap_len_1 + barrier_gate_len + gap_len_2 + plunger_gate_len + gap_len_3 , 0, domain_width, barrier_gate_len/2, 0};
// Boolean fragments
BooleanFragments{ Surface{1:6}; Delete; }{ }

// Extrude upwards to form the gate oxide
// Extrude function has an output that gives all the bottom surfaces
Extrude {0, 0, gate_oxide_thick} {Surface{2:20}; Layers {30};}

// // Physical surfaces for the gate boundaries
Physical Surface("barrier_gate_1_bnd") = {34, 48, 58};
Physical Surface("plunger_gate_1_bnd") = {52, 65, 78};
Physical Surface("barrier_gate_2_bnd") = {72, 81, 87};

// Physical volumes for the gate oxide
Physical Volume("gate_oxide_dot") = {3, 6:19};
Physical Volume("gate_oxide") = {1,2,4,5};

// // Extrude downwards to form the film
Extrude {0, 0, -film_thick} {Surface{2:20}; Layers {30};}

// // Physical surfaces for the source and drain boundaries
Physical Surface("source_bnd") = {91};

// // Physical volumes for the source, drain, and channel
Physical Volume("source") = {20};
Physical Volume("channel") = {24};
Physical Volume("oxide") = {21,23};
Physical Volume("channel_dot") = {26, 28, 31, 34, 36};
Physical Volume("oxide_dot") = {22,25,27,30,33,29,32,35,37,38};

// // Extrude downwards to form the buried oxide
Extrude {0, 0, -box_thick} {
    Surface{97}; Surface{92}; Surface{106}; Surface{106};
    Surface{108}; Surface{115}; Surface{125}; Surface{112}; 
    Surface{122}; Surface{135}; Surface{119}; Surface{132}; 
    Surface{145}; Surface{129}; Surface{142}; Surface{151}; 
    Surface{139}; Surface{148}; Surface{154}; Surface{101}; Layers{30}; 
}
// // Physical surface for the back gate
Physical Surface("back_gate_bnd") = {169,164,160,221,171,178,175,183,186,189,193,196,199,203,206,209,213,216,219};

// // Physical volume for the buried oxide
Physical Volume("buried_oxide") = {39,40,41,42};
Physical Volume("buried_oxide_dot") = {43,46,49,52,55,54,51,48,45,57,44,47,50,53,56};

// // Mesh
Mesh 1;
Mesh 2;
Mesh 3;

//Save
Save "qdfdsoi.msh";
Save "qdfdsoi.xao";