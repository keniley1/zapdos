// Gmsh project created on Thursday Aug 29 15:16:01 2019
dom0Mult = 1e3;
//dom0Mult = 1.0;

lc = 0.001 * dom0Mult;
Point(1111) = {1.068e-5, 0.00101179, 0, 1};

// Point definitions defined in meters
// Lower left
//Point(1) = {0, 0, 0, lc/10.0};

// Upper left
//Point(2) = {0, 0.003, 0, lc/50.0};

// Upper right of boundary
Point(3) = {0.0015 * dom0Mult, 0.002 * dom0Mult, 0.0, lc/6};

// Lower right of boundary (right side of liquid interface)
Point(4) = {0.0015 * dom0Mult, 0, 0, lc/100}; 
Point(43) = {0.0015 * dom0Mult, 0.0002 * dom0Mult, 0.0, lc/10};

Point(6) = {0.0003625 * dom0Mult, 0.002 * dom0Mult, 0, lc/50};
Point(8) = {0.0, 0, 0, lc/100};
Point(81) = {0.0005 * dom0Mult, 0, 0, lc/100};

// For rounded tip of electrode
//lcp = lc / 20000.0;
lcp = lc / 40000.0;
Point(9) = {0.00001 * dom0Mult, 0.00101 * dom0Mult, 0, lcp};
//Point(12) = {-0.00001 * dom0Mult, 0.00101 * dom0Mult, 0, lcp};
Point(12) = {0 * dom0Mult, 0.00100294 * dom0Mult, 0, lcp};
Point(13) = {0, 0.00101356060606 * dom0Mult, 0, lcp};
Circle(85) = {9,13,12};

Point(777) = {0, 0.0009 * dom0Mult, 0, lc/50};
Point(888) = {0, 0.0005 * dom0Mult, 0, lc/20};
Point(889) = {0.0005 * dom0Mult, 0.0005 * dom0Mult, 0, lc/10};

// Create the lines defining the boundaries
// Lines are named according to the points they connect
Line(881) = {8,81}; // Bottom right of domain
Line(814) = {81,4};
//Line(34) = {4,3}; // Right of domain
Line(443) = {4,43};
Line(433) = {43,3};
Line(65) = {6,9}; // Pin lower
Line(36) = {3,6};
Line(12777) = {12,777};
Line(777888) = {777,888};
Line(8888) = {888,8};

//Line(121) = {12,1};
//Line Loop(1) = {4082,828,881,814,34,36,65,85,1260,6030,3040};
//Line Loop(1) = {4082,828,881,814,443,433,36,65,85,1260,6030,304030,403040};
Line Loop(1) = {881,814,443,433,36,65,85,12777,777888,8888};
Plane Surface(0) = {1}; // Define the physical surface of the domain encompassed by curve loop 1

Physical Surface("Inside Domain") = {0};
Physical Curve("bottom") = {881,814};
Physical Curve("top") = {36};
Physical Curve("electrode") = {65,85};
Physical Curve("right") = {34};

lcc = lc/20000;
//Point(999) = {0.000011, 0.00101, 0, lcc};
//Point(121212) = {-0.000011, 0.00101, 0, lcc};
//Point(912) = {0, 0.0010018, 0, lcc};
//Point{999} In Surface{0}; 
//Point{121212} In Surface{0}; 
//Point{912} In Surface{0}; 
//Point(777) = {0, 0.0009 * dom0Mult, 0, lc/100};
//Point(888) = {0, 0.0005 * dom0Mult, 0, lc/20};
//Point(889) = {0.0005 * dom0Mult, 0.0005 * dom0Mult, 0, lc/10};
//Point{777} In Surface{0};
//Point{888} In Surface{0};
Point{889} In Surface{0};
//Point{887} In Surface{0};
Geometry.PointNumbers=1;
Geometry.LineNumbers=1;
Geometry.SurfaceNumbers=1;
//Mesh.ElementOrder=2;
