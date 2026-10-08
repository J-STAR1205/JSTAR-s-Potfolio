// Double Quantum dot in Fully Depleted Silicon On Insulator  (DQDFDSOI)
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
gap_len_4 = 5;    // Length of the gap between barrier 2 and dot 2
gap_len_5 = 5;    // Length of the gap between dot 2 and barrier 3
gap_len_6 = 5;    // Length of the gap between barrier 3 and drain
plunger_gate_len = 15;    // Plunger gate length
barrier_gate_len = 10;     // Barrier gate length
source_drain_len = 20;     // Source/drain length

// z axis dimensions
gate_oxide_thick = 2;
film_thick = 10;
box_thick = 10;

// Rectangle for the simulation domain
channel_len = gap_len_1 + gap_len_2 + gap_len_3 + gap_len_4 + gap_len_5 + gap_len_6 + 3*barrier_gate_len + 2*plunger_gate_len;
domain_len = 2*source_drain_len + channel_len;
Rectangle(1) = {-domain_width/2, -domain_len/2, 0, domain_width, domain_len, 0}; //simulation domain
Printf("%f",domain_len);
// Rectangles for the source, channel and drain
Rectangle(2) = {-channel_width/2, -domain_len/2, 0, channel_width, source_drain_len, 0}; // source area
Rectangle(3) = {-channel_width/2, domain_len/2-source_drain_len, 0, channel_width, source_drain_len, 0}; // drain area
Rectangle(4) = {-channel_width/2, -channel_len/2, 0, channel_width, channel_len, 0}; // channel area

// Rectangles for the plunger gates and barrier gates
gap_list = {gap_len_1, gap_len_2, gap_len_3, gap_len_4, gap_len_5, gap_len_6};
length_list = {plunger_gate_len, barrier_gate_len};
temp = -channel_len/2 + gap_list[0];
For i In {5:9}
    Rectangle(i) = {-domain_width/2, temp  , 0, domain_width, length_list[i%2], 0};
    temp = temp + gap_list[i-4] + length_list[i%2]; 
EndFor

// Boolean fragments
BooleanFragments{ Surface{1:9}; Delete; }{ }

// Extrude upwards to form the gate oxide
// Extrude function has an output that gives all the bottom surfaces
Extrude {0, 0, gate_oxide_thick} {Surface{2:36}; Layers {2};}

// Physical surfaces for the gate boundaries
Physical Surface("barrier_gate_1_bnd") = {55, 69, 79};
Physical Surface("plunger_gate_1_bnd") = {73, 86, 99};
Physical Surface("barrier_gate_2_bnd") = {93, 106, 119};
Physical Surface("plunger_gate_2_bnd") = {113, 126, 139};
Physical Surface("barrier_gate_3_bnd") = {133, 146, 154};

// Physical volumes for the gate oxide
Physical Volume("gate_oxide_dot") = {4, 7:29, 31, 32, 34};
Physical Volume("gate_oxide") = {1, 3, 5, 6, 2, 30, 33, 35};

// Extrude downwards to form the film
Extrude {0, 0, -film_thick} {Surface{2:36}; Layers {2};}

// Physical surfaces for the source and drain boundaries
Physical Surface("source_bnd") = {161};
Physical Surface("drain_bnd") = {163};

// Physical volumes for the source, drain, and channel
Physical Volume("source") = {36};
Physical Volume("drain") = {37};
Physical Volume("channel") = {41, 68};
Physical Volume("oxide") = {38, 40, 65, 70};
Physical Volume("channel_dot") = {43, 45, 48, 51, 54, 57, 60, 63, 66};
Physical Volume("oxide_dot") = {39, 42, 44, 46, 47, 49, 50, 52, 53, 55, 56, 58, 59, 61, 62, 64, 67, 69};

// Extrude downwards to form the buried oxide
Extrude {0, 0, -box_thick} {
    Surface{162}; Surface{172}; Surface{181}; Surface{183};
    Surface{176}; Surface{190}; Surface{200}; Surface{187}; 
    Surface{197}; Surface{210}; Surface{194}; Surface{207}; 
    Surface{220}; Surface{204}; Surface{217}; Surface{230}; 
    Surface{214}; Surface{227}; Surface{240}; Surface{224}; 
    Surface{237}; Surface{234}; Surface{247}; Surface{260}; 
    Surface{244}; Surface{257}; Surface{270}; Surface{254};
    Surface{267}; Surface{275}; Surface{264}; Surface{272};
    Surface{278}; Surface{167}; Surface{250}; Layers{2}; 
}
// Physical surface for the back gate
Physical Surface("back_gate_bnd") = {388, 397, 395, 391, 383, 377, 380,
    383, 373, 370, 367, 356, 359, 363, 399, 352, 349, 339, 342, 345, 335, 332, 329,
    319, 322, 325, 315, 312, 309, 299, 302, 305, 293, 283, 295, 288
};

// Physical volume for the buried oxide
Physical Volume("buried_oxide") = {101, 102, 103, 104, 71, 72, 73, 74};
Physical Volume("buried_oxide_dot") = {75:100, 105};

// Mesh
Mesh 1;
Mesh 2;
Mesh 3;

//Save
Save "dqdfdsoi.msh";
Save "dqdfdsoi.xao";