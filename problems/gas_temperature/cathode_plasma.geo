lc = 1e-2;
lg = 1e-3;
lw = 1e-6;

// Cathode
Point(1) = {0, 0, 0, 1e-3};
Point(2) = {lc, 0, 0, 2e-9};
Line(12) = {1,2};

// Gas
//Point(3) = {, 0, 0, 2e-9};
Point(3) = {lc + lg/2, 0, 0, 20e-6};
Point(4) = {lc + lg, 0, 0, 2e-9};
Line(23) = {2,3};
Line(34) = {3,4};


// Cathode
Physical Line(1) = {12};

// Gas
Physical Line(0) = {23,34};

