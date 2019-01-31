function [r,v] = rv_from_r0v0_ta(r0, v0, dt)
% This function computes the state vector (r,v) from the
% initial state vector (r0,v0) and the change in true anomaly.

% r0 - initial position vector (km)
% v0 - initial velocity vector (km/s)
% dt - change in true anomaly (degrees)

User M-functions required: f_and_g_ta

muo = 398600;

% Compute the f and g functions and their derivatives

[f, g] = f_and_g_ta(r0, v0, dt, muo);
[fdot, gdot] = fDot_and_gDot_ta(r0, v0, dt, muo);
% Compute the final position and velocity vectors:
r = f*r0 + g*v0;
v = fdot*r0 + gdot*v0;
end