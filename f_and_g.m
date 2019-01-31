function [f, g] = f_and_g(x, t, ro, a)

%This function calculates the Lagrange f and g coefficients.

% a: reciprocal of the semimajor axis (1/km)
% ro: the radial position at time to (km)
% t: the time elapsed since ro (s)
% x: the universal anomaly after time t (km^0.5)

muo = 398600;
z = a*x^2;
f = 1 - x^2/ro*stumpC(z);
g = t - 1/sqrt(muo)*x^3*stumpS(z);


function s = stumpS(z)
if (z > 0)
  s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
elseif (z < 0)
  s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
else
  s = 1/6;
end
end

function c = stumpC(z)
if z > 0
  c = (1 - cos(sqrt(z)))/z;
elseif (z < 0)
  c = (cosh(sqrt(-z)) - 1)/(-z);
else
  c = 1/2;
end
end

end