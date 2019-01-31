function [] = GT_from_COE(theta0, a, e, Omega_C0, Omega0, incl, plot_time)

% This code plots the satellite's ground track throughout a user-defined
% period of time. Its inputs are the satellite's orbital elements as well 
% as the period of time over which the ground track is to be plotted.

% theta0: initial true anomaly given by the user in degrees
% a: semi-major axis given by the user in kilometers
% e: orbit's eccentricity
% Omega_C0: initial right ascension of assending node given by the user in degrees
% Omega0: initial argument of perigee given by the user in degrees
% incl: orbit's inclination
% plot_time: the period of time over which the ground track is to be plotted
%%
% Adjusting input data:

theta0 = theta0*pi/180; %in radians
Omega_C0 = Omega_C0*pi/180; %in radians
Omega0 = Omega0*pi/180; %in radians
incl = incl*pi/180; %in radians

% Setting constants:

muo = 398600;
j2 = 0.0010836;
R_Earth = 6378;
w_Earth = 360*(1+1/365.26)/24/3600;
dt = 60; %in seconds
%% 
% Calculating specific angular momentum, orbital period, and the rates of
% change of Omega_C and Omega

h = sqrt(muo*a*(1-e^2));
T = 2*pi/sqrt(muo)*a^(3/2);
Omega_C_dot = -((3*sqrt(muo)*j2*R_Earth^2)/(2*(1-e^2)*a^(7/2)))*cos(incl);
Omega_dot = Omega_C_dot*(5*(sin(incl))^2/2-2)/cos(incl);
%%
% Calculating the eccentric and mean anomaly and the initial time

E0 = 2*atan(tan(theta0/2)*sqrt((1-e)/(1+e)));
M0 = E0-e*sin(E0);
t0 = M0/2/pi*T;

% Preparing necessary empty vectors and setting intials

t = zeros(1,ceil(plot_time/dt)); 
M = zeros(1,length(t)); 
E = zeros(1,length(t));
theta = zeros(1,length(t)); 
Omega_C = zeros(1,length(t)); 
Omega = zeros(1,length(t));
t(1) = t0; 
M(1) = M0; 
E(1) = E0;
theta(1) = theta0;
Omega_C(1) = Omega_C0; 
Omega(1) = Omega0;

coe = {[]}; %classical orbital elements
coe{1} = [h,e,Omega_C(1),incl,Omega(1),theta(1)];

r = {[]}; %position vector in the geocentric equatorial frame
v = {[]}; %velocity vector in the geocentric equatorial frame
r_dash = {[]}; %position vector in the rotating Earth-fixed frame 

theta_Earth = zeros(1,length(t)); %the angle between the stationary and the rotating x-axis
theta_Earth(1) = w_Earth*t(1);

R_theta_Earth = {[]}; %transformation matrix from the geocentric equatorial frame to the rotating Earth-fixed frame
R_theta_Earth{1} = [1,0,0;0,1,0;0,0,1];
ra = zeros(1,length(t)-1); % right ascension 
dec = zeros(1,length(t)-1); %declination
%%
% Calculating data for future time
for i = 1:length(t)-1
    t(i+1) = t(i)+dt;
    M(i+1) = 2*pi*t(i+1)/T;
    E(i+1) = kepler_E(e,M(i+1)); %iterative method to solve the transcendental equation
    
    theta(i+1) = 2*atan(tan(E(i+1)/2)*sqrt((1+e)/(1-e)));
    
    Omega_C(i+1) = Omega_C(i)+Omega_C_dot*dt;
    Omega(i+1) = Omega(i)+Omega_dot*dt;
    coe{i+1} = [h,e,Omega_C(i+1),incl,Omega(i+1),theta(i+1)];
    
    [r{i},v{i}] = sv_from_coe(coe{i+1},muo);
    
    theta_Earth(i+1) = w_Earth*(t(i+1)-t(1));
    R_theta_Earth{i+1} = [cosd(theta_Earth(i+1)),sind(theta_Earth(i+1)),0;-sind(theta_Earth(i+1)),cosd(theta_Earth(i+1)),0;0,0,1];
    
    r_dash{i} = R_theta_Earth{i+1}*r{i}';
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




