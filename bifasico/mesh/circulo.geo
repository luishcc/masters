// Gmsh project created on Mon May 27 15:15:09 2019

Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 0.1};
Point(3) = {-1, 0, 0, 0.1};
Point(4) = {-0.8, 0, 0, 0.4};
Point(5) = {0.8, 0, 0, 0.4};
Point(6) = {-0.9, 0, 0, 0.1};
Point(7) = {0.9, 0, 0, 0.1};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 2};
Circle(3) = {5, 1, 4};
Circle(4) = {4, 1, 5};
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 6};

Line Loop(1) = {1, 2};
Line Loop(2) = {3, 4};
Line Loop(3) = {5, 6};
Plane Surface(1) = {1, 3};
Plane Surface(2) = {2};
Plane Surface(3) = {3, 2};

//+
Physical Line("dirichlet 0") = {1, 2};
//+
Physical Line(2) = {3, 4, 5, 6};
//+
Physical Surface(3) = {1, 2, 3};

