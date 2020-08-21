// This is a geometry file for a unit square.

// Inputs
squareSide = 1;  // arbitrary units; let's call it meters
gridSize = squareSide / 30;

// All numbering counterclockwise from bottom-left corner
Point(1) = {0, 0, 0, gridSize};
Point(2) = {squareSide, 0, 0, gridSize};
Point(3) = {squareSide, squareSide, 0, gridSize};
Point(4) = {0, squareSide, 0, gridSize};

// Now we need to add lines connecting the points
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

// Name the boundaries
Physical Surface(0) = {6};
Physical Curve("bottom") = {1};
Physical Curve("right") = {2};
Physical Curve("top") = {3};
Physical Curve("left") = {4};
Color Red{Physical Surface{0};}

