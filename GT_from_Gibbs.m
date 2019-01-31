function [] = GT_from_Gibbs(R1, R2, R3, plot_time)

% The satellite's ground track is plotted throughout a user-defined period 
% of time where the inputs provided by the user are three position vectors, 
% and the period of time over which the ground track is to be plotted.

% R1, R2, R3: three coplanar geocentric position vectors
% plot_time: the period of time over which the ground track is to be plotted

%%
% Constant:
muo = 398600;
tol = 1e-4; %tolerance for determining if R1, R2 and R3 are coplanar

% Magnitudes of R1, R2 and R3:
r1 = norm(R1);
r2 = norm(R2);
r3 = norm(R3);

% Cross products among R1, R2 and R3:
c12 = cross(R1,R2);
c23 = cross(R2,R3);
c31 = cross(R3,R1);

% Check that R1, R2 and R3 are coplanar; if not set error flag:
if abs(dot(R1,c23)/r1/norm(c23)) > tol
    fprintf('Error: The three position vectors are not coplanoar.');
    return;
end

% Calculate N
N = r1*c23 + r2*c31 + r3*c12;
% Calculate D
D = c12 + c23 + c31;
% Calculate S
S = R1*(r2 - r3) + R2*(r3 - r1) + R3*(r1 - r2);
% Calculate the velocity vector corresponding to the second position vector
V2 = sqrt(muo/norm(N)/norm(D))*(cross(D,R2)/r2 + S);
%%
% Calculating the classical orbital elements (coe) from the state vector (R2,V2)
coe = coe_from_sv(R2,V2);

% The ground track is found from the classical orbital elements
GT_from_COE(coe(6), coe(7), coe(2), coe(3), coe(5), coe(4), plot_time)

end