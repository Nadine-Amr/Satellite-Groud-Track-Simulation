function [tout, yout] = rk4(ode_function, tspan, y0, h)

% This function uses the Runge-Kutta 4 procedure to integrate
% a system of first-order differential equations dy/dt = f(t,y)

% ode_function: handle for M-function in which the derivatives f are computed
% tspan: the vector [t0 tf] giving the time interval for the solution
% y0: column vector of initial values of the vector y
% h: time step
%%
% Calculating the needed constants for Range Kutta-4
n_stages = 4; %the number of points within a time interval that the derivatives are to be computed
a = [0 1/2 1/2 1]; % coefficients for locating the solution points within each time interval
b = [ 0 0 0
1/2 0 0
0 1/2 0
0 0 1]; % coefficients for computing the derivatives at each interior point
c = [1/6 1/3 1/3 1/6]; % coefficients for the computing solution at the end of the time step

%%
t0 = tspan(1);
tf = tspan(2);

t = t0; % t0: initial time
y = y0; % tf: final time

tout = t; % column vector of times at which y was evaluated
yout = y'; % matrix, each row of which contains the components of y evaluated at the corresponding time in tout

while t < tf
   ti = t; % time at the beginning of a time step
   yi = y; % values of y at the beginning of a time step
   
   % Calculate the time derivatives at the n_stages points within the current interval
   for i = 1:n_stages
     t_inner = ti + a(i)*h; % time within a given time step
     y_inner = yi; % values of y within a given time step
     for j = 1:i-1
       y_inner = y_inner + h*b(i,j)*f(:,j);
     end
     f(:,i) = feval(ode_function, t_inner, y_inner);
   end
   h = min(h, tf-t);
   t = t + h;
   y = yi + h*f*c';
   tout = [tout;t]; % adds t to the bottom of the column vector tout
   yout = [yout;y']; % adds y’ to the bottom of the matrix yout
end
end