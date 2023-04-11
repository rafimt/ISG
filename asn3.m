%Assignment 3

clc;
format long;
format compact;
clear;
close all;
fclose all;

%Ex 1

e1 = [1 0 0]';
e2 = [0 1 0]';

g = [1 0; 0 1];

%Ex 2

x1 = 3782970.10; %[m]
y1 = 902154.92;%[m]
z1 = 5038375.59; %[m]

[p1,l1,r1] =xyz2plr(x1,y1,z1); %[rad]
[x2, y2, z2] = plr2xyz(p1,l1,r1); %[m]
[p2,l2,r2] =xyz2plr(x2,y2,z2); %[rad]

rad2deg = 180/pi;

phi_deg = p1*rad2deg %[deg]
lam_deg = l1*rad2deg %[deg]
r_km = r1/1000 %[km]


x_diff = x1/1000-x2/1000 %[km]
y_diff = y1/1000-y2/1000 %[km]
z_diff = z1/1000-z2/1000 %[km]

phi_diff_deg= p1*rad2deg-p2*rad2deg %[deg]

%arc length along a meridian (constant lambda)
R = 6371; %Earth radius[km]
arc_p1 = p1 * R; %[km] (angle should be in rad)
arc_p2 = p2 * R; %[km] (angle should be in rad)
phi_diff_arc = arc_p1-arc_p2; %[km]

lambda_diff_deg = l1*rad2deg-l2*rad2deg %[deg]

arc_l1 = l1 * R * cos(p1); %[km]
arc_l2 = l2 * R * cos(p2); %[km]
lambda_diff_arc = arc_l1-arc_l1 %[km]

radius_diff = r1/1000-r2/1000 %[km]

%Ex3

[B1, L1, h1, e2, a1, b1, N1] = xyz2BLh(x1, y1, z1, 'GRS80');
lat_B1 = B1*rad2deg; %[deg]
long_L1 = L1*rad2deg; %[deg]
height_h1 = h1/1000; %[km]

[x2_2, y2_2, z2_2] = blh2xyz(B1, L1, h1); %[m]
[B2, L2, h2, e2_2, a2, b2, N2] = xyz2BLh(x2_2, y2_2, z2_2, 'GRS80'); 

lat_B2 = B2*rad2deg; %[deg] 
long_L2 = L2*rad2deg; %[deg]
height_h2 = h2/1000; %[km]


x_diff2 = x1/1000-x2_2/1000 %[km]
y_diff2 = y1/1000-y2_2/1000 %[km]
z_diff2 = z1/1000-z2_2/1000 %[km]

B_diff_deg= lat_B1-lat_B2 %[deg]
L_diff_deg = long_L1-long_L2 %[deg]
B_diff_arc = (a1/1000*ellipticE(B1,e2)) - (a2/1000*(ellipticE(B2, e2_2))) %[km] 
L_diff_arc = (L1*N1/1000*cos(B1)) - (L2*N2/1000*cos(B2)) %[km] 
h_diff = height_h1-height_h2 %[km]


%Ex3, b

[B_bessel, L_bessel, h_bessel, e2_bessel, a_bessel, b_bessel, N_bessel] = xyz2BLh(x1, y1, z1, 'BESSEL');
B_bessel_deg = B_bessel*rad2deg %[deg]
B_bessel_arc = distance(0, 0, B_bessel_deg, 0, [a_bessel, sqrt(e2_bessel)])/1000 %km
B_bessel_arc_approx = B_bessel*sqrt(a_bessel*b_bessel)/1000 %km
L_bessel_deg = L_bessel*rad2deg %[deg]
L_bessel_arc = (L_bessel*pi*(N_bessel+h_bessel)*cos(B_bessel))/1000 %[km]
h_bessel_km = h_bessel/1000 %[km]

[B_wgs, L_wgs, h_wgs, e2_wgs, a_wgs, b_wgs, N_wgs] = xyz2BLh(x1, y1, z1, 'WGS84');
B_wgs_deg = B_wgs*rad2deg %[deg]
B_wgs_arc = distance(0, 0, B_wgs_deg, 0, [a_wgs, sqrt(e2_wgs)])/1000 %km
B_wgs_arc_approx = B_wgs*sqrt(a_wgs*b_wgs)/1000 %km
L_wgs_deg = L_wgs*rad2deg %[deg]
L_wgs_arc = (L_wgs*pi*(N_wgs+h_wgs)*cos(B_wgs))/1000 %[km]
h_wgs_km = h_wgs/1000 %[km]

[B_grs, L_grs, h_grs, e2_grs, a_grs, b_grs, N_grs] = xyz2BLh(x1, y1, z1, 'GRS80');
B_grs_deg = B_grs*rad2deg %[deg]
B_grs_arc = distance(0, 0, B_grs_deg, 0, [a_grs, sqrt(e2_grs)])/1000 %km
B_grs_arc_approx = B_grs*sqrt(a_grs*b_grs)/1000 %km
L_grs_deg = L_grs*rad2deg %[deg]
L_grs_arc = (L_grs*pi*(N_grs+h_grs)*cos(B_grs))/1000 %[km]
h_grs_km = h_grs/1000 %[km]

%ratio of difference of bessel and wgs coords compared to grs
ratio_B_bessel = (B_bessel-B_grs)/B_grs*100 %pct ratio
ratio_B_wgs = (B_wgs-B_grs)/B_grs*100 %pct ratio
ratio_L_bessel = (L_bessel-L_grs)/L_grs*100 %pct ratio
ratio_L_wgs = (L_wgs-L_grs)/L_grs*100 %pct ratio
ratio_h_bessel = (h_bessel-h_grs)/h_grs*100 %pct ratio
ratio_h_wgs = (h_wgs-h_grs)/h_grs*100 %pct ratio

%%Exercise 4
diff_B_bessel_deg = B_bessel_deg -phi_deg %deg
diff_B_wgs_deg = B_wgs_deg -phi_deg %deg
diff_B_grs_deg = B_grs_deg -phi_deg %deg

diff_L_bessel_deg = L_bessel_deg - lam_deg %deg
diff_L_wgs_deg = L_wgs_deg - lam_deg %deg
diff_L_grs_deg = L_grs_deg - lam_deg %deg

diff_B_bessel_arc = B_bessel_arc-arc_p1 %km
diff_B_wgs_arc = B_wgs_arc-arc_p1 %km
diff_B_grs_arc = B_grs_arc-arc_p1 %km

diff_L_bessel_arc = L_bessel_arc-arc_l1 %km
diff_L_wgs_arc = L_wgs_arc-arc_l1 %km
diff_L_grs_arc = L_grs_arc-arc_l1 %km

h_sph = r_km-R %km
diff_h_bessel = h_bessel_km - h_sph %km
diff_h_wgs = h_wgs_km - h_sph %km
diff_h_grs = h_grs_km - h_sph %km

dist_3d_bessel = sqrt(diff_B_bessel_arc^2+diff_L_bessel_arc^2+diff_h_bessel^2)
dist_3d_wgs = sqrt(diff_B_wgs_arc^2+diff_L_wgs_arc^2+diff_h_wgs^2)
dist_3d_grs = sqrt(diff_B_grs_arc^2+diff_L_grs_arc^2+diff_h_grs^2)

%%Exercise 5

data= [7205 1492404.605 -4457266.520 4296881.795 -.0155 -.0012 0.0041
7205 1492404.683 -4457266.515 4296881.775 -.0152 -.0014 0.0043
7224 4075539.757 931735.399 4801629.449 -.0160 0.0171 0.0101
7224 4075539.836 931735.313 4801629.400 -.0156 0.0168 0.0104
7232 5085442.772 2668263.635 -2768696.876 -.0015 0.0196 0.0165
7232 5085442.779 2668263.544 -2768696.963 -.0002 0.0193 0.0173
7242 -3950237.046 2522347.621 -4311562.205 -.0393 0.0083 0.0409
7242 -3950236.859 2522347.586 -4311562.417 -.0395 0.0091 0.0415];

X1 = [data(1,2) data(3,2) data(5,2)]';
Y1 = [data(1,3) data(3,3) data(5,3)]';
Z1 = [data(1,4) data(3,4) data(5,4)]';

X2 = [data(2,2) data(4,2) data(6,2)]';
Y2 = [data(2,3) data(4,3) data(6,3)]';
Z2 = [data(2,4) data(4,4) data(6,4)]';

vX1 = [data(1,5) data(3,5) data(5,5)]';
vY1 = [data(1,6) data(3,6) data(5,6)]';
vZ1 = [data(1,7) data(3,7) data(5,7)]';

X1_shift = X1+vX1.*(X2-X1);
Y1_shift = Y1+vY1.*(Y2-Y1);
Z1_shift = Z1+vZ1.*(Z2-Z1);

dX = X2-X1_shift;
dY = Y2-Y1_shift;
dZ = Z2-Z1_shift;

l = [dX' dY' dZ']';

A = [1 0 0 X1_shift(1) 0 Z1_shift(1) -Y1_shift(1)
0 1 0 Y1_shift(1) -Z1_shift(1) 0 X1_shift(1)
0 0 1 Z1_shift(1) Y1_shift(1) -X1_shift(1) 0
1 0 0 X1_shift(2) 0 Z1_shift(2) -Y1_shift(2)
0 1 0 Y1_shift(2) -Z1_shift(2) 0 X1_shift(2)
0 0 1 Z1_shift(2) Y1_shift(2) -X1_shift(2) 0
1 0 0 X1_shift(3) 0 Z1_shift(3) -Y1_shift(3)
0 1 0 Y1_shift(3) -Z1_shift(3) 0 X1_shift(3)
0 0 1 Z1_shift(3) Y1_shift(3) -X1_shift(3) 0];

theta = inv(A'*A)*A'*l;

T = theta(1:3);
D = theta(4);
R = [0 -theta(7) theta(6)
theta(7) 0 -theta(5)
-theta(6) theta(5) 0];

X1_7242 = (data(7,2:4))';
X2_7242 =  (X1_7242+T+D*X1_7242+R*X1_7242)
given_7242 = data(8,2:4)';
coord_diff = (X2_7242-given_7242)
dist_3d = norm(X2_7242-given_7242)

