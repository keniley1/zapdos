// Gmsh project created on Thursday Aug 29 15:16:01 2019
lc = 0.002;

// Point definitions defined in meters
// Lower left
Point(1) = {0, 0, 0, lc/10.0};

// Upper left
//Point(2) = {0, 0.003, 0, lc/50.0};

// Upper right of boundary
Point(3) = {0.0015, 0.002, 0.0, lc/3};

// Lower right of boundary (right side of liquid interface)
Point(4) = {0.0015, 0, 0, lc/5}; 

//Point(5) = {0, 0.001, 0, lc/150.0}; // middle of needle exit (used for meshing)
Point(6) = {0.0003625, 0.002, 0, lc/100};
//Point(7) = {0.0003625, 0.003, 0, lc/50};
Point(8) = {0.0008, 0, 0, lc/10};

// For rounded tip of electrode
Point(9) = {0.00001, 0.00101, 0, lc/10000};
Point(12) = {0, 0.0010029456206381562, 0, lc/10000.0};
Point(13) = {0, 0.00101356060606, 0, 1.0};
Circle(85) = {9,13,12};

// For rounded corner where needle begins to be shaved


// Create the lines defining the boundaries
// Lines are named according to the points they connect
Line(18) = {1,8}; //Bottom left of domain
Line(84) = {8,4}; // Bottom right of domain
Line(34) = {4,3}; // Right of domain
Line(65) = {6,9}; // Pin lower
//Line(76) = {7,6}; // Pin upper
//Line(37) = {3,7}; // Top
Line(36) = {3,6};


Line(121) = {12,1};
//Curve Loop(1) = {34,65,76,37,84,18,121,85};
Curve Loop(1) = {34,65,36,84,18,121,85};
Plane Surface(0) = {1}; // Define the physical surface of the domain encompassed by curve loop 1

Physical Surface("Inside Domain") = {0};
Physical Curve("bottom") = {84,18};
Physical Curve("left_lower") = {121};
Physical Curve("top") = {37};
Physical Curve("right") = {34};
Physical Curve("electrode") = {65,76,85};
//Physical Curve("electrode") = {65,85};
//Physical Curve("top_left") = {76};

Geometry.PointNumbers=1;
Geometry.LineNumbers=1;
Geometry.SurfaceNumbers=1;
//Mesh.ElementOrder=2;
