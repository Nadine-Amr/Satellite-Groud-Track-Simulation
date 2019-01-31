function E = kepler_E(e,M)

% This function uses Newton’s method to solve Kepler’s
% equation E - e*sin(E)= M for the eccentric anomaly,
% given the eccentricity and the mean anomaly.
% e: eccentricity
% M: mean anomaly

% Setting an error tolerance:
error=1.e-8;

% Selecting an initial guess E (the eccentric anomaly in radians):
if M<pi
  E=M+e/2;
else
  E=M-e/2;
end

% Iterate on the non-linear equation of E until E is determined to within the error tolerance:
ratio=1;
while abs(ratio)>error
  ratio=(E-e*sin(E)-M)/(1-e*cos(E));
  E=E-ratio;
end
end