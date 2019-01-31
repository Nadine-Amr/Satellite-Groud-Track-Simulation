function [] = GT_from_R0V0_Lagrange(R0, V0, plot_time)

% This code plots the satellite's ground track throughout a user-defined
% period of time using the Lagrange coefficients. Its inputs are the satellite's initial state vector as well 
% as the period of time over which the ground track is to be plotted.

% R0: initial satellite radius
% V0: initial satellite velocity
% plot_time: the period of time over which the ground track is to be plotted

%%
% Constnts
muo = 398600;
dt = 60; %in seconds
j2 = 0.0010836;
R_Earth = 6378;
w_Earth = 360*(1+1/365.26)/24/3600;
dtheta = 0.5;

% Get the initial classical orbital elements
coe0 = coe_from_sv(R0,V0);
h0 = coe0(1);
e0 = coe0(2);
W0 = coe0(3)*pi/180;
incl0 = coe0(4)*pi/180;
w0 = coe0(5)*pi/180;
theta0 = coe0(6)*pi/180;
a0 = coe0(7);

T = 2*pi/sqrt(muo)*a0^(3/2); % Calculating the orbital period

% Calculating the eccentric and mean anomaly and the initial time
E0 = 2*atan(tan(theta0/2)*sqrt((1-e0)/(1+e0)));
M0 = E0-e0*sin(E0);
t0 = M0/2/pi*T;

% Calculate the transformation matrix from the geocetric equatorial
% coordinates to the perifocal coordinates
QXx = [-sin(W0)*cos(incl0)*sin(w0) + cos(W0)*cos(w0), cos(W0)*cos(incl0)*sin(w0) + sin(W0)*cos(w0), sin(incl0)*sin(w0); -sin(W0)*cos(incl0)*cos(w0) - cos(W0)*sin(w0), cos(W0)*cos(incl0)*cos(w0) - sin(W0)*sin(w0), sin(incl0)*cos(w0); sin(W0)*sin(incl0), -cos(W0)*sin(incl0), cos(incl0)];
QxX = QXx';

rp = {[]}; %position vector in the perifocal frame
vp = {[]}; %velocity vector in the perfocal frame
rp{1} = (QXx*R0')';
vp{1} = (QXx*V0')';

periods = plot_time/T;
theta = zeros(1,ceil(periods*360/dtheta));

theta(1) = t0;

Omega_C = zeros(1,length(theta)); 
Omega = zeros(1,length(theta));
Omega_C(1) = W0; 
Omega(1) = w0;

r = {[]}; %position vector in the geocentric equatorial frame
v = {[]}; %velocity vector in the geocentric equatorial frame
r{1} = R0;
v{1} = V0;
r_dash = {[]}; %position vector in the rotating Earth-fixed frame

theta_Earth = zeros(1,length(theta)); %the angle between the stationary and the rotating x-axis
theta_Earth(1) = w_Earth*theta(1);

R_theta_Earth = {[]}; %transformation matrix from the geocentric equatorial frame to the rotating Earth-fixed frame
R_theta_Earth{1} = [1,0,0;0,1,0;0,0,1];
ra = zeros(1,length(theta)-1); % right ascension 
dec = zeros(1,length(theta)-1); %declination

for i = 1:length(theta)-1
    theta(i+1) = theta(i)+dtheta;
    [rp{i+1},vp{i+1}] = rv_from_r0v0_ta(rp{1},vp{1}, theta(i+1));
    
    r{i+1} = (QxX*rp{i+1}')';
    v{i+1} = (QxX*vp{i+1}')';
    
    coe = coe_from_sv(r{i+1},v{i+1});
    h = coe(1);
    e = coe(2);
    W = coe(3)*pi/180;
    incl = coe(4)*pi/180;
    w = coe(5)*pi/180;
    thetaa = coe(6)*pi/180;
    a = coe(7);
    
    Omega_C_dot = -((3*sqrt(muo)*j2*R_Earth^2)/(2*(1-e^2)*a^(7/2)))*cos(incl);
    Omega_dot = Omega_C_dot*(5*(sin(incl))^2/2-2)/cos(incl);
    
    Omega_C(i+1) = Omega_C(i)+Omega_C_dot*dt;
    Omega(i+1) = Omega(i)+Omega_dot*dt;
    
    W = Omega_C(i+1)*180/pi;
    w = Omega(i+1)*180/pi;
    
    coe = [h e W incl w thetaa a];
    [r{i+1}, v{i+1}] = sv_from_coe(coe,muo);
    
    theta_Earth(i+1) = w_Earth*(theta(i+1)-theta(1));
    R_theta_Earth{i+1} = [cosd(theta_Earth(i+1)),sind(theta_Earth(i+1)),0;-sind(theta_Earth(i+1)),cosd(theta_Earth(i+1)),0;0,0,1];
    
    r_dash{i} = R_theta_Earth{i+1}*r{i+1}';
    [ra(i),dec(i)] = ra_and_dec_from_r(r_dash{i});
end
%%

tol=100; %setting tolerance value to detect when the right ascension has exceeded 360
curve_no=1;
n_curves=1;
k =0;
ra_prev =ra(1);
for i =1:length(ra)
    if abs(ra(i) - ra_prev) > tol
        curve_no=curve_no + 1;
        n_curves=n_curves + 1;
        k =0; %starting a new vector in RA and DEC
    end
    k =k+ 1;
    RA{curve_no}(k) =ra(i);
    Dec{curve_no}(k)=dec(i);
    ra_prev =ra(i);
end
figure(1)
hold on
xlabel('East longitude (degrees)')
ylabel('Latitude (degrees)')
axis equal
grid on
for i =1:n_curves
    plot(RA{i}, Dec{i})
end
axis ([0 360 -90 90])
text( ra(1), dec(1),'Start') %indicating start position
text(ra(end), dec(end),'Finish') %indicating finish position
line([0 360],[0 0],'Color','k') %drawing the equator
title('Ground Track')

figure (2)
plot_3D(r,v)
end