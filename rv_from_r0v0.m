function [R,V] = rv_from_r0v0(R0, V0, t)

% This function computes the state vector (R,V) from the
% initial state vector (R0,V0) and the elapsed time.

% R0: initial position vector (km)
% V0: initial velocity vector (km/s)
% t: elapsed time (s)
%%
muo = 398600;

% Magnitudes of R0 and V0
r0 = norm(R0);
v0 = norm(V0);

% Initial radial velocity
vr0 = dot(R0, V0)/r0;
% Reciprocal of the semimajor axis
alpha = 2/r0 - v0^2/muo;

% The universal anomaly
x = kepler_U(t, r0, vr0, alpha);

% The f and g functions
[f, g] = f_and_g(x, t, r0, alpha);

% The final position vector
R = f*R0 + g*V0;

% The magnitude of R
r = norm(R);

% The derivatives of f and g
[fdot, gdot] = fDot_and_gDot(x, r, r0, alpha);
% the final velocity
V = fdot*R0 + gdot*V0;
end