// Gmsh project created on Thursday Aug 29 15:16:01 2019
d0m = 2e-4; // Domain 0 typical mesh size
d1m = 10e-4; // Domain 1 typical mesh size 
d3m = 0.5e-4; // Domain 3 (interior needle) typical mesh size

needle_mult = 1e0;
gas_mult = 1e1;
dielectric_mult = 4;

// Scaling parameters for spatial domain
d0scale = 1.0;
d1scale = 1.0;

domain_x = 0.03 * d0scale;

dielectric_length = 0.01 * d0scale;
overshoot = 2e-3 * d0scale;  // This is the amount by which the needle is 
                             // longer than the dielectric

gap_distance = 1e-3 * d0scale;

water_depth = 10000e-6 * d0scale;
gas_height = dielectric_length + gap_distance + overshoot;
//needle_thickness = 5e-6 * d0scale;
needle_thickness = 0.25e-3 * d0scale;

domain_y = gas_height + water_depth; 

// Dielectric tube radius and thickness
inner_radius = 4e-3 * d0scale / 2;
outer_radius = 6.35e-3 * d0scale / 2; 
dielectric_thickness = 1.35e-3 * d0scale;

// Water domain
// Right wall is electrode
Point(01) = {0, 0, 0, d1m}; 
Point(02) = {domain_x, 0, 0, d1m};
Point(03) = {domain_x, water_depth, 0, d1m/2}; 
Point(04) = {0, water_depth, 0, d0m/4};

// Gas domain
Point(11) = {needle_thickness, domain_y, 0, d0m * needle_mult};
//Point(12) = {needle_thickness, water_depth + gap_distance, 0, d0m * needle_mult};
//Point(18) = {0, water_depth + gap_distance, 0, d0m * needle_mult};
Point(13) = {domain_x, domain_y, 0, d0m * gas_mult};

// Needle geometry
tip_radius = 1e-4 * d0scale;
Point(112) = {needle_thickness, water_depth + gap_distance + 0.75e-3, 0, d0m};
//Point(113) = {tip_radius, water_depth + gap_distance + tip_radius, d0m};
//Point(114) = {0, water_depth + gap_distance + tip_radius, d0m};
//Circle(115) = {18,114,113};
//Point(113) = {1e-4, 0.0021, 0, 1e-4};
//Point(113) = {11/146000, 297/146000, 0, 1e-4};

//Point(113) = {9.41083e-5, 0.00206618, 0, 1e-4};
//Point(114) = {0, 0.0021, 0, 1e-4};
//Point(18) = {0, 0.002, 0, 1e-4};
//Circle(1115) = {113,114,18};
//Point(113) = {9.41083e-5, 0.00206618, 0, 1e-4};
Point(113) = {9.41083e-05,0.0110662, 0, d0m};
Point(114) = {0, water_depth + gap_distance + tip_radius, 0, d0m};
Point(18) = {0, water_depth + gap_distance, 0, d0m/4};
Circle(1115) = {113,114,18};

// Inner needle geometry point -- upper left
Point(200) = {0, domain_y, 0, d3m};

// Dielectric tubing
Point(14) = {outer_radius, domain_y, 0, d0m * dielectric_mult};
//Point(15) = {outer_radius, water_depth + gap_distance + overshoot, 0, d0m}; 
//Point(16) = {inner_radius, water_depth + gap_distance + overshoot, 0, d0m}; 
// Dielectric tip is a semicircle
radius = (outer_radius - inner_radius) /2;
Point(111) = {inner_radius + (outer_radius - inner_radius)/2, water_depth + gap_distance + overshoot + radius, 0, d0m};
Point(15) = {outer_radius, water_depth + gap_distance + overshoot + radius, 0, d0m * dielectric_mult};
Point(16) = {inner_radius, water_depth + gap_distance + overshoot + radius, 0, d0m * dielectric_mult/2};
Circle(112) = {16,111,15};
Point(17) = {inner_radius, domain_y, 0, d0m * dielectric_mult/2};


// Water domain boundary
Line(12) = {1,2};
Line(23) = {2,3};
Line(34) = {3,4};
Line(41) = {4,1};

// Gas domain boundary
Line(313) = {3,13};
Line(1314) = {13,14};
Line(1415) = {14,15};
Line(1617) = {16,17};
Line(1711) = {17,11};
Line(112113) = {112,113};
Line(11112) = {11,112};
//Line(1112) = {11,12};
//Line(1218) = {12,18};
Line(184) = {18,4};

// Dielectric
Line(1417) = {14,17};

// Needle interior boundary
Line(11200) = {11,200};
Line(20018) = {200,18};


//Line Loop(1) = {-34,313,1314,1415,-112,1617,1711,1112,1218,184};
Line Loop(1) = {-34,313,1314,1415,-112,1617,1711,11112,112113,1115,184};
Plane Surface(0) = {1};

Line Loop(2) = {12,23,34,41};
Plane Surface(1) = {2};

Line Loop(3) = {-1617,112,-1415,1417};
Plane Surface(2) = {3};

Line Loop(4) = {-1115,-112113,-11112,11200,20018};
Plane Surface(3) = {4};

//Physical Curve("electrode") = {1112,1218};
Physical Curve("electrode_side") = {11112};
Physical Curve("electrode_needle") = {112113};
Physical Curve("electrode_tip") = {1115};
Physical Curve("dielectric_left") = {1617};
Physical Curve("dielectric_right") = {1415};
Physical Curve("dielectric_tip") = {112};
Physical Curve("ground") = {23};
Physical Curve("water_bulk") = {12};
Physical Curve("needle_top") = {11200};
Physical Curve("inlet") = {1711};
Physical Curve("water_surface") = {34};
Physical Curve("right") = {313};
Physical Curve("top_right") = {1314};

Physical Surface(0) = {0};
Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};

Geometry.PointNumbers=1;
Geometry.LineNumbers=1;
Geometry.SurfaceNumbers=1;
//+
