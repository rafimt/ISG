clc;
format long;
format compact;
clear;
close all;
fclose all;

%% Assignment 4
%% Exercise 1

%precession 
x_0 = [1 0 0]'; %[m]

yyyy = 1990:2030;
mm = 1; dd = 1; ut1 = 12; minute = 0; second = 0;
jd = zeros(size(yyyy));
delta = zeros(size(yyyy));
alpha = zeros(size(yyyy));
r = zeros(size(yyyy));
for i_yyyy = 1:length(yyyy)
    jd(i_yyyy) = gre2jd(yyyy(i_yyyy), mm, dd, ut1, minute, second);
    P = prec(jd(i_yyyy));
    x_P = P*x_0;
    [delta(i_yyyy), alpha(i_yyyy), r(i_yyyy)] = xyz2plr(x_P(1), x_P(2), x_P(3));
end

alpha_hr = (rad2deg(alpha))/15;
delta_deg = rad2deg(delta);

figure 
[hAx,hLine1,hLine2] = plotyy(yyyy,alpha_hr,yyyy,delta_deg);

title('Precession of ficticious celestial object at (a=0;d=0)')
xlabel('Time (year)')

ylabel(hAx(1),'alpha (h)') % left y-axis 
ylabel(hAx(2),'delta (deg)') % right y-axis

%% Exercise 2
%nutation

mm_n = 1:12; 
jd_n = zeros(length(yyyy),length(mm_n));
for i_yy = 1:length(yyyy)
    for i_mm = 1:length(mm_n)
        jd_n(i_yy,i_mm) = gre2jd(yyyy(i_yy), mm_n(i_mm), dd, ut1, minute, second);
    end
end
jd_nutation = reshape((jd_n'),1,[]);
delta_n = zeros(size(jd_nutation));
alpha_n = zeros(size(jd_nutation));
r_n = zeros(size(jd_nutation));
for i_n = 1:length(jd_nutation)
    N_n = nut(jd_nutation(i_n));
    x_N_n = N_n*x_0;
    [delta_n(i_n), alpha_n(i_n), r_n(i_n)] = xyz2plr(x_N_n(1), x_N_n(2), x_N_n(3));
end

alpha_hr_n = (rad2deg(alpha_n))/15;
delta_deg_n = rad2deg(delta_n);


figure 
[nAx,nLine1,nLine2] = plotyy(jd_nutation,alpha_hr_n,jd_nutation,delta_deg_n);

title('Nutation of ficticious celestial object at (a=0;d=0)')
xlabel('Time (jd)')

ylabel(nAx(1),'alpha (h)') % left y-axis 
ylabel(nAx(2),'delta (deg)') % right y-axis


%% Exercise 3
%21.10.28
xpol_1 = 0.175038/3600; %deg 
ypol_1 = 0.254163/3600; %deg
W_1 = pm(xpol_1,ypol_1);

%21.10.29
xpol_2 = 0.172942/3600; %deg 
ypol_2 = 0.253728/3600; %deg
W_2 = pm(xpol_2,ypol_2);

day = 0:1; 
xpol = [xpol_1 xpol_2]; 
ypol = [ypol_1 ypol_2];
minute = 0:(1/(24*60)):1; 
xpol_interp = interp1(day,xpol,minute);
ypol_interp = interp1(day,ypol,minute);
figure
plot(xpol_interp,ypol_interp)
title('Polar motion interpolated over one full day')
xlabel('xpol')
ylabel('ypol')


W(i_pol) = pm(xpol_1, xpol_2);

%% Exercise 4
lambda = dms2dd(13, 24, 0); %deg
phi = dms2dd(52, 36, 0); %deg
time_ind = ref3d(1)*rot3d(rad2deg(pi/2)-phi,2)*rot3d(lambda,3);

alpha1 = dms2dd(5*15, 8, 42.36351222); %deg 
delta1 = dms2dd(84, 32, 4.5441733);%deg
alpha2 = dms2dd(11*15, 3, 52.22168463); %deg 
delta2 = dms2dd(-53,-57,-0.6966389);%deg
alpha3 = dms2dd(11*15, 13, 58.69508613); %deg 
delta3 = dms2dd(14, 42, 26.9526507);%deg
alpha4 = dms2dd(17, 39, 27.39049431); %deg 
delta4 = dms2dd(49, 55, 3.3683385);%deg
alpha1_r = deg2rad(alpha1); %rad
delta1_r = deg2rad(delta1); %rad
alpha2_r = deg2rad(alpha2); %rad
delta2_r = deg2rad(delta2); %rad
alpha3_r = deg2rad(alpha3); %rad
delta3_r = deg2rad(delta3); %rad
alpha4_r = deg2rad(alpha4); %rad
delta4_r = deg2rad(delta4); %rad
[xi_1(1),xi_1(2),xi_1(3)] = sph2cart(delta1_r,alpha1_r,1); %m
[xi_2(1),xi_2(2),xi_2(3)] = sph2cart(delta2_r,alpha2_r,1); %m
[xi_3(1),xi_3(2),xi_3(3)] = sph2cart(delta3_r,alpha3_r,1); %m
[xi_4(1),xi_4(2),xi_4(3)] = sph2cart(delta4_r,alpha4_r,1); %m

%xpol and ypol interpolated over one day
xpol1 = 0.174245/3600; %deg 
ypol1 = 0.435720/3600; %deg
xpol2 = 0.175463/3600; %deg   
ypol2 = 0.435402/3600; %deg
day = 0:1; 
minute = 0:(1/(24*60)):1-(1/(24*60)); 
xpol = interp1(day,[xpol1 xpol2],minute); %m
ypol = interp1(day,[ypol1 ypol2],minute); %m

%create julian day vector at minute intervals
yr = 2021; mo = 6; dy = 11; hr = 0:23; min = 0:59; sec = 0;
jd = zeros(length(hr),length(min));
for i_hr = 1:length(hr)
    for i_min = 1:length(min)
        jd(i_hr,i_min) = gre2jd(yr,mo,dy,hr(i_hr),min(i_min),sec);
    end
end
jd = reshape((jd'),1,[]); %jd

p1 = zeros(size(jd));
l1 = zeros(size(jd));
for i_n = 1:length(jd)
    P = prec(jd(i_n));
    N = nut(jd(i_n));
    GMST = gmst(jd(i_n)); %hours
    W = pm(xpol(i_n),ypol(i_n));
    xg = time_ind*W*rot3d(GMST*15,3)*N*P*xi_1'; %m
    [p1(i_n), l1(i_n),~] = cart2sph(xg(1), xg(2), xg(3));
    l1(i_n) = 90-rad2deg(l1(i_n)); %deg
    p1(i_n) = rad2deg(p1(i_n)); %deg
    if p1(i_n) < 0
        p1(i_n) = p1(i_n)+360;
    end
end

p2 = zeros(size(jd));
l2 = zeros(size(jd));
for i_n = 1:length(jd)
    P = prec(jd(i_n));
    N = nut(jd(i_n));
    GMST = gmst(jd(i_n));
    W = pm(xpol(i_n),ypol(i_n));
    xg = time_ind*W*rot3d(GMST*15,3)*N*P*xi_2';
    [p2(i_n), l2(i_n),~] = cart2sph(xg(1), xg(2), xg(3));
    l2(i_n) = 90-rad2deg(l2(i_n));
    p2(i_n) = rad2deg(p2(i_n));
    if p2(i_n) < 0
        p2(i_n) = p2(i_n)+360;
    end
end

p3 = zeros(size(jd));
l3 = zeros(size(jd));
for i_n = 1:length(jd)
    P = prec(jd(i_n));
    N = nut(jd(i_n));
    GMST = gmst(jd(i_n));
    W = pm(xpol(i_n),ypol(i_n));
    xg = time_ind*W*rot3d(GMST*15,3)*N*P*xi_3';
    [p3(i_n), l3(i_n),~] = cart2sph(xg(1), xg(2), xg(3));
     l3(i_n) = 90-rad2deg(l3(i_n));
     p3(i_n) = rad2deg(p3(i_n));
    if p3(i_n) < 0
        p3(i_n) = p3(i_n)+360;
    end
end

p4 = zeros(size(jd));
l4 = zeros(size(jd));
for i_n = 1:length(jd)
    P = prec(jd(i_n));
    N = nut(jd(i_n));
    GMST = gmst(jd(i_n));
    W = pm(xpol(i_n),ypol(i_n));
    xg = time_ind*W*rot3d(GMST*15,3)*N*P*xi_4';
    [p4(i_n), l4(i_n),~] = cart2sph(xg(1), xg(2), xg(3));
     l4(i_n) = 90-rad2deg(l4(i_n));
     p4(i_n) = rad2deg(p4(i_n));
    if p4(i_n) < 0
        p4(i_n) = p4(i_n)+360;
    end
end

hours = 0:(24/1440):24-(24/1440);

figure 
subplot(2,2,1)
[nAx,nLine1,nLine2] = plotyy(hours,p1,hours,l1);
title('Source 0454+844')
xlabel('Time (hours)')
xticks([0 6 12 18 24])
xticklabels({'0','6','12','18','24'})
ylabel(nAx(1),'Azimuth (deg)') % left y-axis 
ylabel(nAx(2),'Elevation (deg)') % right y-axis
subplot(2,2,2)
[nAx,nLine1,nLine2] = plotyy(hours,p2,hours,l2);
title('Source 1101?536')
xlabel('Time (hours)')
xticks([0 6 12 18 24])
xticklabels({'0','6','12','18','24'})
ylabel(nAx(1),'Azimuth (deg)') % left y-axis 
ylabel(nAx(2),'Elevation (deg)') % right y-axis
subplot(2,2,3)
[nAx,nLine1,nLine2] = plotyy(hours,p3,hours,l3);
title('Source: 1111+149')
xlabel('Time (hours)')
xticks([0 6 12 18 24])
xticklabels({'0','6','12','18','24'})
ylabel(nAx(1),'Azimuth (deg)') % left y-axis 
ylabel(nAx(2),'Elevation (deg)') % right y-axis
subplot(2,2,4)
[nAx,nLine1,nLine2] = plotyy(hours,p4,hours,l4);
title('Source: 1738+499')
xlabel('Time (hours)')
xticks([0 6 12 18 24])
xticklabels({'0','6','12','18','24'})
ylabel(nAx(1),'Azimuth (deg)') % left y-axis 
ylabel(nAx(2),'Elevation (deg)') % right y-axis


figure

hsky1 = skyplot(p1,l1, 'm')
hold on
hsky2 = skyplot(p2,l2, 'g')
hold on
hsky3 = skyplot(p3,l3,'r')
hold on
hsky4 = skyplot(p4,l4,'b')
h = zeros(4, 1);
h(1) = plot(NaN,NaN,'om');
h(2) = plot(NaN,NaN,'og');
h(3) = plot(NaN,NaN,'or');
h(4) = plot(NaN,NaN,'ob');
legend(h, '0454+844','1101?536','1111+149','1738+499')
title('Celestial objects in local horizon system')
