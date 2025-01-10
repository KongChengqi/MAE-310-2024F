%% Matlab mesh
%% quarter-plate-with-hole-quad_1, Created by Gmsh
%% ASCII
clear msh;
msh.nbNod = 32;
msh.POS = [
1 -1 0;
1 1 0;
-1 1 0;
-0.7 -1 0;
-1 -0.7 0;
-0.7878679656440357 -0.7878679656440357 0;
-1 -0.4571428571432273 0;
-1 -0.2142857142867692 0;
-1 0.02857142856968109 0;
-1 0.2714285714267745 0;
-1 0.5142857142846611 0;
-1 0.7571428571423739 0;
-0.7500000000006948 1 0;
-0.5000000000013898 1 0;
-0.250000000002084 1 0;
-2.756350703236876e-12 1 0;
0.2499999999979157 1 0;
0.4999999999986104 1 0;
0.7499999999993052 1 0;
1 0.7500000000006948 0;
1 0.5000000000013898 0;
1 0.250000000002084 0;
1 2.756350703236876e-12 0;
1 -0.2499999999979157 0;
1 -0.4999999999986104 0;
1 -0.7499999999993052 0;
0.7571428571432381 -1 0;
0.5142857142867491 -1 0;
0.271428571430319 -1 0;
0.02857142857322559 -1 0;
-0.214285714284661 -1 0;
-0.4571428571423739 -1 0;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES =[
 4 6 12
 5 7 11
 7 8 11
 8 9 11
 9 10 11
 10 11 11
 11 12 11
 12 3 11
 3 13 9
 13 14 9
 14 15 9
 15 16 9
 16 17 9
 17 18 9
 18 19 9
 19 2 9
 2 20 8
 20 21 8
 21 22 8
 22 23 8
 23 24 8
 24 25 8
 25 26 8
 26 1 8
 1 27 10
 27 28 10
 28 29 10
 29 30 10
 30 31 10
 31 32 10
 32 4 10
];
