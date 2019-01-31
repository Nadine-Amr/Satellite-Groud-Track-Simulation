function dfdt = rates(t,f)

% This function calculates first and second time derivatives of the orbit's radius
% Dr: velocity (r')
% D2r: acceleration (r")
% f: column vector containing r and Dr at time t
% dfdt: column vector containing Dr and D2r at time t
%%
% Constant
muo = 398600;

%%
r = [f(1);f(2);f(3)];
% Calculating norm of radius
rnorm = norm(r);
Dr = [f(4);f(5);f(6)];

for i=1:3
    D2r(i) = -muo*r(i)/(rnorm^3);
end
D2r = D2r';

%%
dfdt = [Dr; D2r];
end