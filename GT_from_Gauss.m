function [] = GT_from_Gauss(azim, elev, LST, lat, R1, R2, R3, t1, t2, t3, plot_time)
% The satellite's ground track is plotted throughout a user-defined period 
% of time where the inputs provided by the user are three ground track locations, 
% the angles of observation, and the period of time over which the ground track 
% is to be plotted. The Earth's oblateness and rotation are taken into consideration.

% azim: vector of the azimuths at each of the ground stations
% elev: vector of the elevations at each of the ground stations
% LST: vector of the local sideral times at each of the ground stations
% lat: vector of the latitudes at each of the ground stations
% R1, R2, R3: the observation site position vectors at t1, t2, t3
% plot_time: the period of time over which the ground track is to be plotted

%%
% Constants
muo = 398600;

for i=1:3
  % Calculating direction cosines in topocentric horizon coordinates
  l = [sind(azim(i))*cosd(elev(i)); cosd(azim(i))*sind(elev(i)); sind(elev(i))];
  % Calculating rotation matrix to topocentric equatorial coordinates
  QxX = [-sind(LST(i)), -cosd(LST(i))*sind(lat(i)), cosd(LST(i))*cosd(lat(i)); cosd(LST(i)), -sind(LST(i))*sind(lat(i)), sind(LST(i))*cosd(lat(i)); 0, cosd(lat(i)), sind(lat(i))];
  % Calculating direction cosines in topocentric equatorial coordinates
  Rho(i) = QxX*l;
end

Rho1 = Rho(1);
Rho2 = Rho(2);
Rho3 = Rho(3);

% Find time intervals between observations
tau1 = t1 - t2;
tau3 = t3 - t2;
tau = tau3 - tau1;

% Obtain cross products among the direction cosine vectors
p1 = cross(Rho2,Rho3);
p2 = cross(Rho1,Rho3);
p3 = cross(Rho1,Rho2);

Do = dot(Rho1,p1); %scalar triple product of Rho1, Rho2 and Rho3
D = [[dot(R1,p1) dot(R1,p2) dot(R1,p3)]
[dot(R2,p1) dot(R2,p2) dot(R2,p3)]
[dot(R3,p1) dot(R3,p2) dot(R3,p3)]]; %matrix of the nine scalar triple products of R1, R2 and R3 with p1, p2 and p3

E = dot(R2,Rho2);

%constants in the expression relating slant range to geocentric radius
A = 1/Do*(-D(1,2)*tau3/tau + D(2,2) + D(3,2)*tau1/tau);
B = 1/6/Do*(D(1,2)*(tau3^2 - tau^2)*tau3/tau + D(3,2)*(tau^2 - tau1^2)*tau1/tau);
%coefficients of the 8th order polynomial in the estimated geocentric radius x
a = -(A^2 + 2*A*E + norm(R2)^2);
b = -2*muo*B*(A + E);
c = -(muo*B)^2;

%Calculate the roots of the 8th order polynomial rho1, rho2, rho3 - the slant ranges at t1, t2, t3
% r1, r2, r3 - the position vectors at t1, t2, t3 using MATLAB’s polynomial 'roots' solver:
Roots = roots([1 0 a 0 0 b 0 0 c]);
%Find the positive real root:
x = posroot(Roots);

%Initial estimates for Lagrange coefficients
f1 = 1 - 1/2*muo*tau1^2/x^3;
f3 = 1 - 1/2*muo*tau3^2/x^3;
g1 = tau1 - 1/6*muo*(tau1/x)^3;
g3 = tau3 - 1/6*muo*(tau3/x)^3;

%the slant ranges at t1, t2, t3
rho2 = A + muo*B/x^3;
rho1 = 1/Do*((6*(D(3,1)*tau1/tau3 + D(2,1)*tau/tau3)*x^3 + muo*D(3,1)*(tau^2 - tau1^2)*tau1/tau3) /(6*x^3 + muo*(tau^2 - tau3^2)) - D(1,1));
rho3 = 1/Do*((6*(D(1,3)*tau3/tau1 - D(2,3)*tau/tau1)*x^3 + muo*D(1,3)*(tau^2 - tau3^2)*tau3/tau1) /(6*x^3 + muo*(tau^2 - tau1^2)) - D(3,3));

%Initial estimates
r1 = R1 + rho1*Rho1;
r2 = R2 + rho2*Rho2;
r3 = R3 + rho3*Rho3;
v2 = (-f3*r1 + f1*r3)/(f1*g3 - f3*g1);

%We iterate to obtain better estimates
rho1_old = rho1; 
rho2_old = rho2; 
rho3_old = rho3;
diff1 = 1; 
diff2 = 1; 
diff3 = 1;
n = 0;
nmax = 1000; %max number of iterations
tol = 1.e-8; %the error tolerance determining convergence

while ((diff1 > tol) && (diff2 > tol) && (diff3 > tol) && (n < nmax))
  n = n+1;
%Compute quantities required by universal Kepler’s equation
  ro = norm(r2); %magnitude of the position vector
  vo = norm(v2); %magnitude of the velocity vector
  vro = dot(v2,r2)/ro; %radial velocity component
  a = 2/ro - vo^2/muo; %reciprocal of the semimajor axis
%Solve universal Kepler’s equation at times tau1 and tau3 for universal anomalies x1 and x3:
  x1 = kepler_U(tau1, ro, vro, a);
  x3 = kepler_U(tau3, ro, vro, a);
%Calculate the Lagrange f and g coefficients at times tau1 and tau3:
  [ff1, gg1] = f_and_g(x1, tau1, ro, a);
  [ff3, gg3] = f_and_g(x3, tau3, ro, a);
%Update the f and g functions at times tau1 and tau3 by averaging old and new:
  f1 = (f1 + ff1)/2;
  f3 = (f3 + ff3)/2;
  g1 = (g1 + gg1)/2;
  g3 = (g3 + gg3)/2;
 
  c1 = g3/(f1*g3 - f3*g1);
  c3 = -g1/(f1*g3 - f3*g1);
  rho1 = 1/Do*( -D(1,1) + 1/c1*D(2,1) - c3/c1*D(3,1));
  rho2 = 1/Do*( -c1*D(1,2) + D(2,2) - c3*D(3,2));
  rho3 = 1/Do*(-c1/c3*D(1,3) + 1/c3*D(2,3) - D(3,3));
  
  r1 = R1 + rho1*Rho1;
  r2 = R2 + rho2*Rho2;
  r3 = R3 + rho3*Rho3;
  
  v2 = (-f3*r1 + f1*r3)/(f1*g3 - f3*g1);

  diff1 = abs(rho1 - rho1_old);
  diff2 = abs(rho2 - rho2_old);
  diff3 = abs(rho3 - rho3_old);

%Update the slant ranges:
  rho1_old = rho1; 
  rho2_old = rho2; 
  rho3_old = rho3;
end

%State vector for the central observation:
r = r2;
v = v2;
%%
% Calculating the classical orbital elements (coe) from the state vector (r,v)
coe = coe_from_sv(r,v);

% The ground track is found from the classical orbital elements
GT_from_COE(coe(6), coe(7), coe(2), coe(3), coe(5), coe(4), plot_time)

end