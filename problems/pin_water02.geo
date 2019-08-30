// Gmsh project created on Thursday Aug 29 15:16:01 2019
lc = 0.0005;

// Point definitions defined in meters
// Lower left
Point(1) = {0, 0, 0, lc/200.0};

// Upper left
//Point(2) = {0, 0, 0.003, lc/50.0};

// Upper right of boundary
Point(3) = {0, 0.005, 0.003, lc};

// Lower right of boundary (right side of liquid interface)
Point(4) = {0, 0.005, 0, lc/5}; 

//Point(5) = {0, 0, 0.001, lc/150.0}; // middle of needle exit (used for meshing)
Point(6) = {0, 0.0003625, 0.002, lc/10};
Point(7) = {0, 0.0003625, 0.003, lc/10};
Point(8) = {0, 0.0008, 0, lc/40};

// For rounded tip of electrode
Point(5) = {0, 0.00001, 0.00101, lc/600};
Point(9) = {0, 0, 0.001, lc/600.0}; // middle of needle exit (used for meshing)
Point(10) = {0, 0.000005, 0.001002, lc/600.0};
// Create the lines defining the boundaries
// Lines are named according to the points they connect
//Line(15) = {5,1}; // Left, below capillary tube
Line(19) = {9,1};
//Line(52) = {2,5}; // Left, above capillary tube
//Line(14) = {1,4}; // Bottom of domain (liquid interface)
Line(18) = {1,8}; //Bottom left of domain
Line(84) = {8,4}; // Bottom right of domain
Line(34) = {4,3}; // Right of domain
//Line(32) = {3,2}; // Top
Line(65) = {6,5}; // Pin lower
Line(76) = {7,6}; // Pin upper
Line(37) = {3,7}; // Top

BSpline(85) = {5, 10, 9};

// connect the loops encompassing the domain
//Curve Loop(1) = {15,52,14,34,32};
//Curve Loop(1) = {15,14,34,65,76,37};
Curve Loop(1) = {34,65,76,37,84,18,19,85};
Plane Surface(1) = {1}; // Define the physical surface of the domain encompassed by curve loop 1

Physical Surface("Inside Domain") = {1};
//Physical Curve("bottom") = {14};
//Physical Curve("bottom_left") = {18};
//Physical Curve("bottom_right") = {84};
Physical Curve("bottom") = {84,18};
//Physical Curve("left_lower") = {15};
Physical Curve("left_lower") = {19};
//Physical Curve("left_upper") = {52};
//Physical Curve("top") = {32};
Physical Curve("top") = {37};
Physical Curve("right") = {34};
//Physical Curve("electrode") = {65,76};
Physical Curve("electrode") = {65,76,85};

Point(20) = {0, 0.000015, 0.00101, lc/1000};
Point{20} In Surface{1}; 
Geometry.PointNumbers=1;
Geometry.LineNumbers=1;
Geometry.SurfaceNumbers=1;
//+
//Circle(85) = {5, 10, 9};
//+
//BSpline(85) = {5, 10, 9};
//Spline(85) = {5,10,9};

