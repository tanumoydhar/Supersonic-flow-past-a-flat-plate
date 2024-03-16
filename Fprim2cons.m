
function [F1, F2, F3, F5] = Fprim2cons(rho, u, p, v, T, mu, lambda, k, c_v, dx, dy, call_case)
% Compute the flux vector F from the primitive flow field variables

% Preallocate other arrays with zeros 
[numy,numx] = size(u);
F1 = zeros(numy,numx);
F2 = zeros(numy,numx);
F3 = zeros(numy,numx);
F5 = zeros(numy,numx);


F1 = rho.*v; % From continuity
tau_yx = xyshear(u, v, mu, dx, dy, call_case);
    

F2 = rho.*u.*v - tau_yx; % From x-momentum


tau_yy = yyshear(u, v, lambda, mu, dx, dy, call_case);
F3 = rho.*v.^2 + p - tau_yy; % From y-momentum



q_y = heat_flux_y(T, k, dy, call_case);
F5 = (rho.*(c_v*T + (u.^2 + v.^2)/2) + p).*v - u.*tau_yx - v.*tau_yy + q_y;