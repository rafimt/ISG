%% Assignment 2

clc;
format long g;
format compact;
clear;
close all;
fclose all;

%% Exercise 1
% To compute a 3D rotation matrix, we require two inputs: the rotation angle 
% and rotation axis. The output of this function is a 3x3 rotation matrix. 
% 
% This function results in the identity matrix when rotated by 2pi (or 360 
% degrees). 
% 
% 360 degrees are used here because this function takes degrees as an input. 

axis = 2;
Rotation_2pi = rot3d(360, axis)
%% 
% Rotating by 0 would similarly result in the identity matrix, as would 
% rotating by pi*pi. 

Rotation_0 = rot3d(0,axis)
Rotation_pi_times_pi = rot3d(180,axis)*rot3d(180,axis)
%% 
% For a left-handed system, positive angles rotate the vector in the clockwise 
% direction.
% 
% Thus, the first rotation of -90 degrees over the z-axis results in base 
% vector e2. This could also be achieved by rotating over the same axis by +270 
% degrees. 
% 
%  Then, a rotation of -90 degrees over the x-axis results in base vector 
% e3. 

x1 = [1, 0, 0]'
x2_neg = rot3d(-90,3)*x1
x2_pos = rot3d(270,3)*x1
x3 = rot3d(-90,1)*x2_neg
%% Exercise 2
% Squaring a 3D reflection matrix results in the identity matrix. Multiplying 
% the 3D reflection matrices going around the x, y, and z axes results in a negative 
% identity matrix. 

Refl_squared = ref3d(1)*ref3d(1)
Refl_times_refl = ref3d(1)*ref3d(2)*ref3d(3)
%% Exercise 3
% The equation for computing the eccentric anomaly E from a given mean anomaly 
% M and eccentricity e is overdetermined, thus requires an iterative solution. 
% 
% The input values for the function that iterates to find E include M and 
% e, but it could also be extended to include the stop criterion rather than specifying 
% it within the function itself.
% 
%  


%grid

%M from 0 to 2pi, e from 0 to 1)
M_an = 0:0.01:2*pi;
ecc_1 = 0:0.001:(1-0.001);

%compute E-m at each grid point

%initialize matrices of ones to speed up computation
E_an = ones(length(M_an),length(ecc_1));
E_M_diff = ones(length(M_an),length(ecc_1));
i_ter = ones(length(M_an),length(ecc_1));

for i_M_an = 1:length(M_an)
    for i_ecc_1 = 1:length(ecc_1)
        [E_an(i_M_an,i_ecc_1),i_ter(i_M_an,i_ecc_1)]=kepler(M_an(i_M_an),ecc_1(i_ecc_1));
        E_M_diff(i_M_an,i_ecc_1)=E_an(i_M_an,i_ecc_1)-M_an(i_M_an);
    end
end

%plot

%% 
% The following figure shows the difference between the eccentric anomaly 
% and mean anomaly. The largest difference occurs around M=1.5 and e=0.999. The 
% differences generally tend to get larger as the eccentricity grows. 

figure
grid on
mesh(E_M_diff)
xlabel('e')
ylabel('M')
zlabel('E-M')
colorbar
%% 
% The following figure shows the iterations required to achieve a solution 
% that fulfills the stop criterion. The most iterations needed occur at the same 
% point as the spike in the above plot showing the greatest difference between 
% E and M. 

figure
grid on
mesh(i_ter)
xlabel('e')
ylabel('M')
zlabel('Iterations')
colorbar
%% Exercise 4
% To convert from decimal degrees into hour, minute, second, the degrees would 
% have to be divided by 15. Otherwise the conversion would remain unchanged. 360 
% degrees is comparable to 24 hours, so each hour is 15 degrees.
% 
% However, going from hour angle decimal to hour minute second would not 
% require any conversion because the hour would already be in the correct unit. 
% We would just have to update to robustness check for the hour unit to go from 
% 0 to 24 rather than -360 to 360. 
% 
% For decimal degrees, the sign should remain negative for the degree, but 
% not the minutes and seconds. This is achieved by taking the negative integer 
% value for the degree as is, then taking the absolute values of the degree to 
% do the calculations for minutes and seconds. 
% 
% Going back the other way, the minutes and seconds are either added or subtracted 
% back to the degree to get the decimal degree, depending on the sign of the degree. 
% 
% The following is an example of the conversion from -22.342 decimal degrees 
% to degree, minute, second, and then back to decimal degrees. 

[d, m, s] = dd2dms(-22.342)
[dd] = dms2dd(d, m, s)
%% 
% 
%% Exercise 5
% In CET, summer time is from the last Sunday in March to the last Sunday in 
% October (UTC+2:00). The given time is in winter time, so it is +1 ahead of UTC. 
% 
% The dUT1 value obtained from the IERS website form 2021-11-22 is -0.1073087, 
% so this is added to the UTC time to convert to UT1. 
% 
% The resulting Julian Day is as follows. 


%The following parameters are for the corresponding UT1
yyyy = 2021;
mm = 11;
dd = 21;
hour = 23;
minute = 59;
second = 59.8926913;

jdUT1 = gre2jd(yyyy, mm, dd, hour, minute, second)
%% Exercise 5a
% The Earth Rotation Angle is found by inserting the Julian Day into a linear 
% equation, and then ensuring that the ERA remains between 0 and 2pi. Finally, 
% it is converted from radians into decimal degrees, then into degree, minute, 
% second. 
% 
% The resulting earth rotation angle is 60 deg, 55', 24''. 

tUT1 = jdUT1 - 2451545.0;
ERA = 2*pi*(0.7790572732640+1.00273781191135448*tUT1); 

%Ensure that ERA is between 0 and 2*pi

while (ERA <= 0 || ERA > 2*pi)
    if ERA <= 0
        ERA = ERA + (2*pi); %rad
    elseif ERA > 2*pi
        ERA = ERA - (2*pi); %rad
    end
end
 
ERA_deg = ERA*(180/pi); %deg
[ERA_d, ERA_m, ERA_s] = dd2dms(ERA_deg); %dms
ERA_s_round = round(ERA_s); %seconds
ERA_d
ERA_m
ERA_s_round
%% Exercise 5b
% GAST is determined by taking the Julian Day, computing GMST, and then adding 
% the Equation of Equinox. 


t = (jdUT1 - 2451545.0)/36525;
GMST = (F(jdUT1)*86400+24110.54841-86400/2+8640184.812866*t+0.093104*t*t-6.2e-6*t*t*t)/3600; %hour

while (GMST <= 0 || GMST > 24)
    if GMST <= 0
        GMST = GMST + 24;
    elseif GMST > 24
        GMST = GMST - 24;
    end
end
GMST
%% 
% The delta psi value -112.77 obtained from IERS is in milliarcseconds, 
% so it must be converted to radians.
% 
% The resulting GAST following this transformation is 04:04:49. 


delta_psi = -112.77*4.8481368e-9; %milliarcsec to rad
eps = (84381.448-46.8150*t-(5.9e-4)*t*t+(1.813e-3)*t*t*t)*(pi/648000); %rad
omega = (450160.28-482890.539*t+7.455*t*t+0.008*t*t*t)*(pi/648000)+2*pi*F(-5*t);%rad
EqE = delta_psi*cos(eps)+((2.64e-3)*sin(omega)+(6.3e-5)*sin(2*omega))*(pi/648000); %rad
EqE = EqE*(12/pi); %hour
GAST = GMST+EqE; %hour

while (GAST <= 0 || GAST > 24)
    if GAST <= 0
        GAST = GAST + 24;
    elseif GAST > 24
        GAST = GAST - 24;
    end
end

[GAST_h, GAST_m, GAST_s] = dd2dms(GAST);
GAST_s = round(GAST_s);
GAST_h
GAST_m
GAST_s

%% Exercise 6
% The sun contributes by far the most out of the individual gravitational terms. 
% It contribute about 99.97% to the total gravitational term.


%convert milliarcsecond to second
dTT = 20e-3; %s

%defining constant LG
LG = 6.969290134e-10;

dTCG = (1+(LG/(1-LG)))*dTT; %s

%Vector containing the gravitational terms
GM = [1.32712440018e20 2.2032e13 3.24859e14 4.9048695e12 4.282837e13 6.26325e10 1.26686534e17 3.7931187e16 5.793939e15 6.836529e15 8.71e11 1.108e12]'; %m3s-2

%Vector containing the mean distances to earth
d = [149597870.7 91691000 41400000 385000.6 78340000 414000000 628730000 1275000000 2723950000 4351400000 5890000000 96.1*149597870.7]'*1000; %m

V_ext = sum(GM./d);
%Percent contribution of the gravitational terms of individual solar system
%bodies
names = {'Sun'; 'Mercury'; 'Venus'; 'Moon'; 'Mars'; 'Ceres'; 'Jupiter'; 'Saturn'; 'Uranus'; 'Neptune'; 'Pluto'; 'Eris';};
grav_pct = ((GM./d)./V_ext)*100; %pct 
[grav_pct_planets] = arrayfun(@(x) [names{x} ' ' num2str(grav_pct(x))], 1:length(names),'un',0)'
%% 
% The following ratio of the gravitational relativistic term over the special 
% relativistic term results in 1.97, which is approximately double. Special relativity 
% takes into account the speed of light, and does not consider the effects of 
% gravity. 

v_e = 30e3; %ms-1
%speed of light 
c = 299792458; %m/s
ratio = V_ext/(v_e^2/2)
%% 
% The equivalent group delay of 20 ms in TT to TCB is still approximately 
% 20 ms. The time difference is very small, in the order of 10^-10 seconds. The 
% baseline distance is about 9.3 cm. 


dTCB = dTCG/(1-(1/c^2)*((v_e^2)/2+V_ext));%s
d_time = abs(dTCB-dTT) %s
d_dist = d_time*c %m