Geometry.AutoCoherence = 0;

R = 0.5;
L = 2.0;

Small = 1.0;

Point(1) = {L, -L, 0,Small};
Point(2) = {L, L, 0,Small};
Point(3) = {-L, L, 0,Small};
Point(4) = {-L, -L, 0,Small};
Point(5) = {-L + R, -L, 0,Small};
Point(6) = {-L, -L + R, 0,Small};
Point(7) = {-L + Cos(Pi/4) * R, -L + Sin(Pi/4) * R, 0,Small};

Circle(1) = {5, 4, 7};
Circle(2) = {7, 4, 6};

Line(3) = {6, 3};
Line(4) = {3, 2};
Line(5) = {2, 1};
Line(6) = {1, 5};


Curve Loop(1) = {4, 5,6,1, 2, 3};
Plane Surface(1) = {1};

Transfinite Surface{1};
Recombine Surface{1};
Mesh.ElementOrder = 1;
Mesh.Algorithm = 8;





//+
Physical Curve("right", 8) = {5};
//+
Physical Curve("top", 9) = {4};



//+
Physical Curve("down", 10) = {6};
//+
Physical Curve("left", 11) = {3};
//+
Physical Curve("circle", 12) = {1};
