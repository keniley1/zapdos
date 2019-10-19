// Gmsh project created on Thursday Aug 29 15:16:01 2019
//dom0Mult = 1e3;
dom0Mult = 1.0;

lc = 0.0025 * dom0Mult;

Point(3) = {0.002 * dom0Mult, 0.0015 * dom0Mult, 0.0, lc/6};

// Lower right of boundary (right side of liquid interface)
Point(4) = {0.002 * dom0Mult, 0, 0, lc/1000}; 
Point(43) = {0.002 * dom0Mult, 0.0002 * dom0Mult, 0.0, lc/10};

Point(6) = {0.0005 * dom0Mult, 0.0015 * dom0Mult, 0, lc/100};
//Point(8) = {0.0, 0, 0, lc/2500000};
Point(8) = {0.0, 0, 0, lc/625000};
Point(81) = {0.0008 * dom0Mult, 0, 0, lc/4000};
Point(82) = {0.8e-5 * dom0Mult, 1e-6, 0, lc/1000};
Point(83) = {0.14e-5 * dom0Mult, 0.6e-6, 0, lc/1000};

// For rounded tip of electrode
//lcp = lc / 20000.0;
lcp = lc / 40000.0;
//Point(9) = {0.00001 * dom0Mult, 0.00101 * dom0Mult, 0, lcp*6};
//Point(12) = {-0.00001 * dom0Mult, 0.00101 * dom0Mult, 0, lcp};
Point(12) = {0 * dom0Mult, 0.001 * dom0Mult, 0, lcp*5};
Point(13) = {0, 0.0015 * dom0Mult, 0, lcp};
Circle(85) = {6,13,12};

Point(777) = {0, 0.0009 * dom0Mult, 0, lc/200};
Point(888) = {0, 0.0005 * dom0Mult, 0, lc/40};
Point(889) = {0.0005 * dom0Mult, 0.0005 * dom0Mult, 0, lc/15};
Point(890) = {0.0012 * dom0Mult, 0.00055 * dom0Mult, 0, lc/5};

// Go up!
//Point(661) = {0.0003625 * dom0Mult, 0.0025 * dom0Mult, 0, lc/100};
//Point(663) = {0.002 * dom0Mult, 0.0025 * dom0Mult, 0, lc/10};


// Create line for meshing 
Point(800) = {0 * dom0Mult, 1e-5 * dom0Mult, 0.0, lc / 400};
Point(400) = {0.002 * dom0Mult, 1e-5 * dom0Mult, 0.0, lc / 200};
Line(800400) = {800,400};

// Create the lines defining the boundaries
// Lines are named according to the points they connect
Line(881) = {8,81}; // Bottom right of domain
Line(814) = {81,4};
//Line(34) = {4,3}; // Right of domain
//Line(443) = {4,43};
Line(4400) = {4,400};
Line(40043) = {400,43};
Line(433) = {43,3};
//Line(65) = {6,9}; // Pin lower
//Line(3663) = {3,663};
//Line(663661) = {663,661};
//Line(6616) = {661,6};
Line(36) = {3,6};
Line(12777) = {12,777};
Line(777888) = {777,888};
//Line(8888) = {888,8};
Line(888800) = {888,800};
Line(8008) = {800,8};



// Water Boundary
Point(40) = {0.002 * dom0Mult, -1.5e-5 * dom0Mult, 0, lc/400};
Point(80) = {0.0, -1.5e-5 * dom0Mult, 0, lc/800};
Line(440) = {4,40};
Line(4080) = {40,80};
Line(808) = {80,8};

//Line Loop(5) = {881,814,443,433,3663,663661,6616,65,85,12777,777888,8888};
//Line Loop(5) = {881,814,4400,40043,433,3663,663661,6616,65,85,12777,777888,888800,8008};
Line Loop(5) = {881,814,4400,40043,433,36,85,12777,777888,888800,8008};
Plane Surface(0) = {5}; // Define the physical surface of the domain encompassed by curve loop 1

Line Loop(1) = {881,814,440,4080,808};
//Line Loop(1) = {440,4080,808};
Plane Surface(1) = {1};

Physical Surface(0) = {0};
//Physical Curve("interface") = {881,814};
//Physical Curve("top") = {663661};
Physical Curve("top") = {36};
//Physical Curve("electrode") = {6616,65,85};
Physical Curve("electrode") = {85};
Physical Curve("right") = {443,433};

Physical Surface(1) = {1};
Physical Curve("bottom") = {4080};
Physical Curve("bottom_right") = {440};
Physical Curve("bottom_left") = {808};

//Extrude { Surface{1}; Layers{ {1,1}, {0.1,0.3}}; }

// Make structured
//Recombine Surface{0};
//Recombine Surface{1};

lcc = lc/20000;
Point{889} In Surface{0};
Point{890} In Surface{0};
Point{82} In Surface{0};
//Point{83} In Surface{0};
Line{800400} In Surface{0};
Geometry.PointNumbers=2;
Geometry.LineNumbers=2;
Geometry.SurfaceNumbers=2;
//Mesh.ElementOrder=2;
