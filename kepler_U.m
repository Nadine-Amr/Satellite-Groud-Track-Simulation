function x = kepler_U(dt, ro, vro, a)

% This function uses Newton's method to solve the universal Kepler equation for the universal anomaly.

% dt: time since x = 0 (s)
% ro: radial position (km) when x = 0
% vro: radial velocity (km/s) when x = 0
% a: reciprocal of the semimajor axis (1/km)

muo = 398600;
%Set error tolerance and a limit on the number of iterations:
error = 1.e-8;
nMax = 10000;

%Initial value for x (the universal anomaly):
x = sqrt(muo)*abs(a)*dt;

%Iterate until convergence occurs 
n = 0;
ratio = 1;
while (abs(ratio) > error && n <= nMax)
  n = n + 1;
  C = stumpC(a*x^2);
  S = stumpS(a*x^2);
  F = ro*vro/sqrt(muo)*x^2*C + (1 - a*ro)*x^3*S + ro*x - sqrt(muo)*dt;
  dFdx = ro*vro/sqrt(muo)*x*(1 - a*x^2*S) + (1 - a*ro)*x^2*C + ro;
  ratio = F/dFdx;

  x = x - ratio;
end

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

