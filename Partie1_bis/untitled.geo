//+
Point(1) = {0.0, 0.0, 0, 1.0};
//+
Point(2) = {0.0, 1.0, 0, 1.0};
//+
Point(3) = {1.0, 1.0, 0, 1.0};
//+
Point(4) = {1.0, 0.0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Physical Curve("right", 5) = {1};
//+
Physical Curve("top", 6) = {2};
//+
Physical Curve("left", 7) = {3};
//+
Physical Curve("bot", 8) = {4};
