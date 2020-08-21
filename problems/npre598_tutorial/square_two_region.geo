// Inputs
squareSide = 1;  // arbitrary units; let's call it meters
meshThickness = squareSide / 10;
gridSize = squareSide / 20;

// All numbering counterclockwise from bottom-left corner
Point(1) = {0, 0, 0, gridSize};
Point(2) = {squareSide, 0, 0, gridSize};
Point(3) = {squareSide, squareSide, 0, gridSize};
Point(4) = {0, squareSide, 0, gridSize};

// Add points at x = 0, dividing the domain into 2 equal subdomains
Point(5) = {squareSide/2, 0, 0, gridSize};
Point(6) = {squareSide/2, squareSide, 0, gridSize};

// Now we need to add lines connecting the points
Line(1) = {1, 5};
Line(2) = {5, 6};
Line(3) = {6, 4};
Line(4) = {4, 1};
Line(5) = {5, 2};
Line(6) = {2, 3};
Line(7) = {3, 6};
Curve Loop(6) = {1, 2, 3, 4};
Curve Loop(7) = {5, 6, 7, -2};
Plane Surface(7) = {6};
Plane Surface(8) = {7};

// Name the boundaries
Physical Surface(0) = {7};
Physical Surface(1) = {8};
Physical Curve("bottom_left") = {1};
Physical Curve("bottom_right") = {5};
Physical Curve("top_left") = {3};
Physical Curve("top_right") = {7};
Physical Curve("left") = {4};
Physical Curve("right") = {6};

// Note that there is no definition of the interfacial boundary at x = 0. 
// This will be defined in the MOOSE input file.

