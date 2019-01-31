function [r,v] = sv_from_coe(coe,muo)

% This function computes the state vector (r,v) from the
% classical orbital elements (coe)

% muo: gravitational parameter (km^3/s^2)
% coe: the classical orbital elements [h e RA incl w TA] where
% h = specific angular momentum (km^2/s)
% e = eccentricity
% RA = right ascension of the ascending node (rad)
% incl = inclination of the orbit (rad)
% w = argument of perigee (rad)
% TA = true anomaly (rad)

%%
h=coe(1);
e=coe(2);
RA=coe(3);
incl=coe(4);
w=coe(5);
TA=coe(6);

% Calculatng the position and velocity vectors in the perifocal frame
rp=(h^2/muo)*(1/(1+e*cos(TA)))*(cos(TA)*[1;0;0]+sin(TA)*[0;1;0]);
vp=(muo/h)*(-sin(TA)*[1;0;0]+(e+cos(TA))*[0;1;0]);

% Calculating the rotation matrices
R3_W=[cos(RA),sin(RA),0;-sin(RA),cos(RA),0;0,0,1]; %rotation matrix about the z-axis through the angle w
R1_i=[1,0,0;0,cos(incl),sin(incl);0,-sin(incl),cos(incl)]; %rotation matrix about the x-axis through the angle incl
R3_w=[cos(w),sin(w),0;-sin(w),cos(w),0;0,0,1]; %rotation matrix about the z-axis through the angle RA
Q_pX=(R3_w*R1_i*R3_W)'; %matrix of the transformation from perifocal to geocentric equatorial frame

% Calculating the poition and velocity vectors in the geocentric equatorial frame
r=Q_pX*rp;
v=Q_pX*vp;
r=r';
v=v';
end