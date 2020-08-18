// Inputs
squareSide = 1;  // arbitrary units; let's call it meters
meshThickness = squareSide / 10;
gridSize = squareSide / 20;

// All numbering counterclockwise from bottom-left corner
Point(1) = {0, 0, 0, gridSize};
Point(2) = {squareSide*2.5, 0, 0, gridSize};
Point(3) = {squareSide*2.5, squareSide, 0, gridSize};
Point(4) = {squareSide, squareSide, 0, gridSize};
Point(5) = {squareSide, squareSide*3, 0, gridSize};
Point(6) = {0, squareSide*3, 0, gridSize};

// Now we need to add lines connecting the points
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Curve Loop(5) = {1, 2, 3, 4, 5, 6};
Plane Surface(6) = {5};

// Name the boundaries
Physical Surface(0) = {6};
Physical Curve("bottom") = {1};
Physical Curve("lower_right") = {2};
Physical Curve("lower_top") = {3};
Physical Curve("upper_right") = {4};
Physical Curve("upper_top") = {5};
Physical Curve("left") = {6};
