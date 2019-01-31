function [ra,dec] = ra_and_dec_from_r(r)

% This function calculates the right ascension and the
% declination from the geocentric equatorial position vector.
% r: geocentric equatorial position vector


% Calculating l, m, n (the direction cosines of r)
l = r(1)/norm(r);
m = r(2)/norm(r);
n = r(3)/norm(r);

dec=asind(n); % dec: declination in degrees

if m > 0
   ra=acosd(l/cosd(dec)); % ra: right ascension in degrees
else
   ra=360 - acosd(l/cosd(dec)); % ra: right ascension in degrees
end