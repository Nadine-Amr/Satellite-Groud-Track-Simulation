function [fdot, gdot] = fDot_and_gDot(x, r, ro, a)

% This function calculates the time derivatives of the Lagrange f and g coefficients.

% a: reciprocal of the semimajor axis (1/km)
% ro: the radial position at time to (km)
% r: the radial position after time t (km)
% x: the universal anomaly after time t (km^0.5)
%%
muo = 398600;
z = a*x^2;

fdot = sqrt(muo)/r/ro*(z*stumpS(z) - 1)*x;

gdot = 1 - x^2/r*stumpC(z);
end