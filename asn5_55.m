%% Assignment 5
%% Exercise 1

clc; format compact; clear; close all; fclose all;
fid = fopen('data.txt');
C = textscan(fid, '%s %f %f %f %f %f %f %s');
fclose(fid);

%all vals in deg except name, e, n(day-1), ut1
[name,i,omega,e,w,M1,n,ut1]=C{1:8};

%% 
% The ISS and Sentinel 1B are in low earth orbits, and it takes about 90 
% minutes for them to circle the earth. Molniya 1T is also an LEO, but it is at 
% times quite far from the earth, resulting in an orbital period of about 7.5 
% hours. Galileo9 is in a medium earth orbit, with an orbital period of about 
% 14 hours. Beidou IGSO 3 and Meteosat 11 both have orbital periods of 1 day, 
% as geosynchronous and geostationary satellites respectively. They have a mean 
% motion of 1, meaning that they corotate with the earth. 

%orbital period
T_day = 1./n; %day
T = T_day.*24; %hour decimal
[T_hr, T_min, T_sec] = dd2dms(T); %hour decimal to hms
orbitalperiod = table(name, T_hr, T_min) %rounded to integer min

%% Exercise 2
% The perigee and apogee of all satellites except Molniya are identical to the 
% semi-major axis because the eccentricities are 0. Only for Molniya, with an 
% eccentricity of 0.62, the perigee and apogee vary, meaning that its distance 
% to earth varies depending on where it is in the orbit. 

%%Exercise 2
T_seconds = T.*3600; %seconds
GM = 398600.44; %km3s-2
a = nthroot(((T_seconds.^2.*GM)./(4*pi^2)),3); %[km]
a_meters = fix(a.*1000); %m
perigee = fix(a.*1000.*(1-e)); %[m]
apogee = fix(a.*1000.*(1+e)); %[m]
a_perigee_apogee = table(name, a_meters, perigee, apogee)
%% Exercise 3

E = ones(6,1440);
M = ones(6,1440);
r_dist = ones(6,1440);
v = ones(6,1440);
r_x = ones(6,1440);
r_y = ones(6,1440);
r_z = ones(6,1440);

for sat=1:6
    for i_n = 1:1440 %minutes per day
        t2 = (i_n-1)*1; %[min]
        t_day2 = t2/1440; %[day]
        M(sat,i_n) = deg2rad(M1(sat))+(t_day2*n(sat)*2*pi); %[rad]
        [E(sat,i_n),~] = kepler(M(sat,i_n),e(sat),1e-9); %[rad]  
        r_dist(sat,i_n) = a(sat)*(1-e(sat)*cos(E(sat,i_n))); %km
        v(sat,i_n) = atan2(sqrt(1-e(sat)^2)*sin(E(sat,i_n)),cos(E(sat,i_n))-e(sat)); %[rad]
        r_orb = r_dist(sat,i_n)*[cos(v(sat,i_n)) sin(v(sat,i_n)) 0]'; %km
        r_x(sat,i_n) = r_orb(1); %km
        r_y(sat,i_n) = r_orb(2); %km
        r_z(sat,i_n) = r_orb(3); %km

    end   
end

r_ecix = ones(size(r_x));
r_eciy = ones(size(r_x));
r_eciz = ones(size(r_x));
for sat=1:6
    for mins=1:1440
        r_eci = rot3d(-omega(sat),3)*rot3d(-i(sat),1)*rot3d(-w(sat),3)*[r_x(sat,mins) r_y(sat,mins) r_z(sat,mins)]';
        r_ecix(sat,mins) = r_eci(1);
        r_eciy(sat,mins) = r_eci(2);
        r_eciz(sat,mins) = r_eci(3);
    end
end
%% 
% Plotting the orbits in the ECI system allows us to see the orbital motion 
% in a frame that originates at the center of the earth but is fixed with respect 
% to the stars. It is a useful depiction to visualize parameters such as the shape 
% of the orbit and distance from the earth. The inclination is also clearly visible, 
% for example the Meteosat 11 satellite moves in a very low inclination, close 
% to the equator. The Sentinel 1B on the other hand has an orbit that is very 
% close to the poles.

figure;
clf;
hold on
view(3);
[x,y,z] = sphere(180);
x=6371.*x;
y=6371.*y; 
z=6371.*z;
surf(x,y,z,'FaceColor', 'none','EdgeColor', 0.5*[1 1 1]);
for sat=1:6
pl(sat)=plot3(r_ecix(sat,:),r_eciy(sat,:),r_eciz(sat,:),'DisplayName',char(name(sat)));
end
xlabel("x (km)");
ylabel("y (km)");
zlabel("z (km)");
legend(pl)
title('24 hour ECI orbit');
grid on;
hold off;
v_mag = ones(size(r_x));
v_x = ones(size(r_x));
v_y = ones(size(r_x));
v_z = ones(size(r_x));
%velocity
for mins=1:1440
    for sat=1:6
        p = r_dist(sat,mins)*(1+e(sat)*cos(v(sat,mins)));
        C = sqrt(p*GM);
        C_hat = deg2rad([sind(omega(sat))*sind(i(sat)) -cosd(omega(sat))*sind(i(sat)) cosd(i(sat))]');
        K_hat = deg2rad([cosd(omega(sat)) sind(omega(sat)) 0]');
        P = deg2rad(cosd(w(sat))*K_hat+sind(w(sat))*(cross(C_hat,K_hat)));
        Q = deg2rad(-sind(w(sat))*K_hat+cosd(w(sat))*(cross(C_hat,K_hat)));
        v_orb = (C/p)*(-sin(v(sat,mins))*P+(e(sat)+cos(v(sat,mins)))*Q);
        v_mag(sat,mins) = norm(v_orb);
        v_x(sat,mins) = v_orb(1);%km/s
        v_y(sat,mins) = v_orb(2); %km/s
        v_z(sat,mins) = v_orb(3); %km/s
        
    end
end
    
    
%% 
% The figures below show the velocities of the satellites that oscillate 
% within a period that matches the orbital period. 
% 
% Aside from Molniya, all of the satellites have a constant velocity. The 
% speed at which they travel is the same throughout the orbit because they have 
% a circular orbit. Only Molniya has an inconstant velocity due to its elliptical 
% orbit.
% 
% For Molniya, the velocity is highest at the perigee and lowest at the apogee, 
% because the satellite travels faster when it is closer to the earth. 

for sat=1:6
    
    figure
    hold on
    x_hr = 0:(24/1440):24-(24/1440);
    plot(x_hr,v_x(sat,:))
    plot(x_hr,v_y(sat,:))
    plot(x_hr,v_z(sat,:))
    title(sprintf('24h ECI velocity components of %s',char(name(sat))));
    xlabel('Time [h]');
    xticks([0 6 12 18 24])
    ylabel('Velocity [km/s]');
    legend('Velocty X','Velocity Y', 'Velocity Z');
    figure
    plot(x_hr,v_mag(sat,:))
    title(sprintf('24h ECI velocity magnitude of %s',char(name(sat))));
    xlabel('Time [h]');
    xticks([0 6 12 18 24])
    ylabel('Velocity [km/s]');

end
%% Exercise 4


%%Exercise 4

dt = datetime(string(ut1),'InputFormat','yyyy,MM,dd,HH,mm,ss');    
r_ecefx = ones(size(r_x));
r_ecefy = ones(size(r_x));
r_ecefz = ones(size(r_x));
jd = zeros(24,60,6);

for sat = 1:6
hours = dt.Hour(sat)+((1:24)-1)*1; 
minutes = dt.Minute(sat)+((1:60)-1)*1;
    for i_hr = 1:length(hours)
        for i_min = 1:length(minutes)
            jd(i_hr,i_min,sat)= gre2jd(dt.Year(sat),dt.Month(sat),dt.Day(sat),hours(i_hr),minutes(i_min),dt.Second(sat));
        end
    end
end
jd_min = reshape(permute(jd,[2,1,3]), [],6);
jd_min=jd_min';
for sat = 1:6
    for i_n = 1:1440
        GMST = gmst(jd_min(sat,i_n));
        r_ecef = rot3d(GMST*15,3)*[r_ecix(sat,i_n) r_eciy(sat,i_n) r_eciz(sat,i_n)]';
        r_ecefx(sat,i_n) = r_ecef(1);
        r_ecefy(sat,i_n) = r_ecef(2);
        r_ecefz(sat,i_n) = r_ecef(3);
    
    end
end
%% 
% The ECEF system similarly has its origin at the center of the earth, but 
% it rotates with the earth with respect to the stars. Thus, we can see, for example, 
% the orbits of the low earth satellites moving around the earth many times in 
% one day, or the Beidou satellite's fixed latitudinal orbit. Some satellites 
% appear to have global coverage, while others are not visible from some parts 
% of the earth within the course of one day.

figure;
clf;
view(3);
hold on;
[x,y,z] = sphere(180); x=6371.*x; y=6371.*y; z=6371.*z;
earth=surf(x,y,-z,'FaceColor','none','EdgeColor',[1 1 1]);
cdata=imread('earth.jpeg');
set(earth,'FaceColor','texturemap','CData',cdata,'FaceAlpha',1,'EdgeColor','none');
axis equal; axis auto; axis vis3d;
for sat=1:6
pl2(sat)=plot3(r_ecefx(sat,:),r_ecefy(sat,:),r_ecefz(sat,:),'DisplayName',char(name(sat)));
end
xlabel("x (km)"); ylabel("y (km)"); zlabel("z (km)");
legend(pl2);
title('24 hour ECEF orbit');
grid on; hold off;

%% Exercise 5

%%Exercise 5

phi=ones(size(r_x));
lam=ones(size(r_x));
for sat=1:6
    for i_n = 1:1440
        [phi(sat,i_n),lam(sat,i_n),~]=xyz2plr(r_ecefx(sat,i_n).*1000,r_ecefy(sat,i_n).*1000,r_ecefz(sat,i_n).*1000);
        lam(sat,i_n) = rad2deg(lam(sat,i_n));
        phi(sat,i_n) = rad2deg(phi(sat,i_n));
    end
end
%% 
% The satellite ground tracks provide another perspective of the satellites 
% in the ECEF system. We can see Beidou's orbit in a figure-8 shape over China, 
% Monliya's orbit that speeds up when it is nearer to the earth over the Northern 
% hemisphere, Meteosat 11 is stationary off the coast of Africa. It is also interesting 
% to compare the ISS and Sentinel 1B, the latter having much better coverage in 
% the polar regions. 

figure;
clf;
cdata=imread('earth.jpeg');
imagesc([-180 180],[90,-90],cdata);
set(gca,'ydir','normal');
hold on;
for sat=1:6
    pl3(sat)=plot(lam(sat,:),phi(sat,:),'.','DisplayName',char(name(sat)));
end
xlabel("latitude");
ylabel("longitude");
legend(pl3)
title('ECEF satellite ground tracks');
grid on;
hold off;
%% Exercise 6


%%Exercise 6

r_berx = 3782971.01/1000; %km
r_bery = 902152.49/1000; %km
r_berz = 5038375.59/1000; %km
R_K = [r_berx r_bery r_berz]; %km
[phi_ber,lam_ber,~] = xyz2plr(r_berx,r_bery,r_berz); %rad
d = ones(6,1440);
for sat =1:6
    for i_n = 1:1440
        d(sat,i_n) = (magnitude(R_K,[r_ecefx(sat,i_n) r_ecefy(sat,i_n) r_ecefz(sat,i_n)]))*1000; %m
    end
end
%% 
% The minimum distance of 465150 meters occurred with the ISS at 18:53. 

min_dist=min(d(:));
min_dist_m = fix(min_dist) %m
[row_dist,col_dist] = find(d==min_dist);
min_dist_sat = name(row_dist)
time_hour = fix(col_dist/60)
time_min = fix(60*((col_dist/60)-time_hour))
%% 
% The minimum zenith distance of 92 arc minutes occurred also with the ISS 
% at the same time.

zenith = ones(size(r_x));
for sat=1:6
    for i_n = 1:1440
        r_topo = rad2deg(M(sat,i_n))*rot3d(90-rad2deg(phi_ber),2)*rot3d(rad2deg(lam_ber),3)*[r_ecefx(sat,i_n) r_ecefy(sat,i_n) r_ecefz(sat,i_n)]';
        zenith(sat,i_n) = 90-(atand(r_topo(3)/(sqrt(r_topo(1)^2+r_topo(2)^2)))); %deg
    end
end

min_zenith=min(zenith(:));
min_zenith_arcmin = fix(min_zenith*60)
[row_zenith,col_zenith] = find(zenith==min_zenith);
min_zenith_sat = name(row_zenith)
time_hr_zenith = fix(col_zenith/60)
time_min_zenith = fix(60*((col_zenith/60)-time_hr_zenith))
%% 
% This makes sense given that the ISS is very close to the earth, and the 
% tracks have fairly large global coverage (aside from the polar regions). The 
% Sentinel 1B has similar orbital parameters but a more polar orbit, and it did 
% not approach Berlin as much as the ISS did. The satellite ground tracks plotted 
% above show the ISS flying very close overhead to Berlin. Beidou and Molniya 
% do not go near Berlin because they are intended to orbit more over China and 
% Russia respectively. The Galielo satellite did not approach Berlin at this particular 
% epoch. The Meteosat 11 is very far away and geostationary over the equator, 
% so it would not approach Berlin.