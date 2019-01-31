function [] = GT_from_TLE(filename, plot_time)

% This code plots the satellite's ground track throughout a user-defined
% period of time. Its inputs are the name of the satellite's two-line 
% element file as well as the period of time over which the ground track is 
% to be plotted.

% filename: name of the satellite's two-line element file
% plot_time: the period of time over which the ground track is to be plotted
%%
% Constant
muo = 398600;

% Open the TLE file specified by the user
fid = fopen(filename, 'r');

% Extract the two lines
L1 = fscanf(fid,'%c',69);
garbage = fscanf(fid, '%c', 2);
L2 = fscanf(fid,'%c',69);
fclose(fid);

% Extract direct orbital elements
incl = str2double(L2(9:16)); % inclination in degrees
RA = str2double(L2(18:25)); % right ascension of ascending node in degrees
e = str2double(L2(27:33))/1e7; % eccentricity
w = str2double(L2(35:42)); % argument of perigee in degrees

% Calsulate indirect orbital elements
n = str2double(L2(53:63)); % mean motion in revolutions per day
a = (muo/(n*2*pi/(24*3600))^2)^(1/3); % semi-major axis in kilometers
h = sqrt(a*muo*(1-(e^2))); % specific angular momentum in km^2/s

M = str2double(L2(44:51));  % mean anomaly in degrees
E = kepler_E(e, (M*pi/180)); % eccentric anomaly in radians
TA = 2*atan(sqrt((1+e)/(1-e))*tan(E/2)); % true anomaly in radians

TA = TA*180/pi; % true anomaly in degrees

coe = [h e RA incl w TA a];

%%
% The ground track is found from the classical orbital elements
GT_from_COE(coe(6), coe(7), coe(2), coe(3), coe(5), coe(4), plot_time)
end