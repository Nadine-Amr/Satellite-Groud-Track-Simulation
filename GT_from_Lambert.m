function [] = GT_from_Lambert(R1, R2, t, string, plot_time)

% The satellite's ground track is plotted throughout a user-defined period 
% of time where the inputs provided by the user are two position vectors, 
% the time interval between them, and the period of time over which the 
% ground track is to be plotted.

% R1, R2: initial and final position vectors
% t: the time between R1 to R2 in seconds
% string: this is 'pro' if the orbit is prograde & 'retro' if the orbit is retrograde
% plot_time: the period of time over which the ground track is to be plotted
%%
% Constant:
muo = 398600;

% Magnitudes of R1 and R2:
r1 = norm(R1);
r2 = norm(R2);

c12 = cross(R1, R2); % cross product of R1 into R2
theta = acos(dot(R1,R2)/r1/r2); %angle between R1 and R2

% The orbit is assumed to be prograde if the user did not correctly define this
if (nargin < 4 || (~strcmp(string,'retro') && (~strcmp(string,'pro'))))
    string = 'pro';
end

% The value of theta is corrected
if strcmp(string,'pro')
  if c12(3) <= 0
    theta = 2*pi - theta;
  end
elseif strcmp(string,'retro')
  if c12(3) >= 0
    theta = 2*pi - theta;
  end
end

A = sin(theta)*sqrt(r1*r2/(1 - cos(theta))); %a defined constant

% We determine approximately where F(z,t) changes sign to get an initial condition for z
z = -100;
while F(z,t) < 0
  z = z + 0.1;
end

% Setting an error tolerance and a limit on the number of iterations
tol = 1.e-8;
nmax = 50000;

% Iterate on Equation 5.45 until z is determined to within the error tolerance:
ratio = 1;
n = 0;
while ((abs(ratio) > tol) && (n <= nmax))
  n = n + 1;
  ratio = F(z,t)/dFdz(z);
  z = z - ratio;
end

% Calculate Lagrange coefficients
f = 1 - y(z)/r1;
g = A*sqrt(y(z)/muo);

% Calculate the velocity corresponding to R1
V1 = 1/g*(R2 - f*R1);
%%
% Calculating the classical orbital elements (coe) from the state vector (R1,V1)
coe = coe_from_sv(R1,V1);

% The ground track is found from the classical orbital elements
GT_from_COE(coe(6), coe(7), coe(2), coe(3), coe(5), coe(4), plot_time)
%%
% Some of the needed functions
function dum = y(z)
dum = r1 + r2 + A*(z*S(z) - 1)/sqrt(C(z));
end

function dum = F(z,t)
dum = (y(z)/C(z))^1.5*S(z) + A*sqrt(y(z)) - sqrt(muo)*t;
end

function dum = dFdz(z)
if z == 0
  dum = sqrt(2)/40*y(0)^1.5 + A/8*(sqrt(y(0)) + A*sqrt(1/2/y(0)));
else
  dum = (y(z)/C(z))^1.5*(1/2/z*(C(z) - 3*S(z)/2/C(z)) + 3*S(z)^2/4/C(z)) + A/8*(3*S(z)/C(z)*sqrt(y(z)) + A*sqrt(C(z)/y(z)));
end
end

function dum = C(z)
dum = stumpC(z);
end

function dum = S(z)
dum = stumpS(z);
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
