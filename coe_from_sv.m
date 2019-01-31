function coe = coe_from_sv(R,V)

% This function computes the classical orbital elements (coe)
% from the state vector (R,V)

% R: position vector in the geocentric equatorial frame (km)
% V: velocity vector in the geocentric equatorial frame (km)
%%

%Constants
muo = 398600;
eps = 1.e-10; %value below which the eccentricity is considered to be zero

r = norm(R);
v = norm(V);
vr = dot(R,V)/r; %radial velocity component 
H = cross(R,V); %angular momentum vector
h = norm(H);

incl = acos(H(3)/h); %calculating inclination

N = cross([0 0 1],H); %calculating the node line vector
n = norm(N);

if n ~= 0
  RA = acos(N(1)/n); %calculating the right ascension of the ascending node in radians
  if N(2) < 0
    RA = 2*pi - RA;
  end
else
  RA = 0;
end

E = 1/muo*((v^2 - muo/r)*R - r*vr*V); %calculating the eccentricity vector
e = norm(E);

if n ~= 0
  if e > eps
    w = acos(dot(N,E)/n/e); %calculating the argument of perigee in radians
    if E(3) < 0
      w = 2*pi - w;
    end
  else
    w = 0;
  end
else
  w = 0;
end

if e > eps
  TA = acos(dot(E,R)/e/r); %calculating the true anomaly in radians
  if vr < 0
    TA = 2*pi - TA;
  end
else
  cp = cross(N,R);
  if cp(3) >= 0
    TA = acos(dot(N,R)/n/r);
  else
    TA = 2*pi - acos(dot(N,R)/n/r);
  end
end

a = h^2/muo/(1 - e^2); %calculating the semimajor axis

incl = incl*180/pi;
RA = RA*180/pi;
w = w*180/pi;
TA = TA*180/pi;

coe = [h e RA incl w TA a]; %constructing a vector of orbital elements [h e RA incl w TA a]
end 