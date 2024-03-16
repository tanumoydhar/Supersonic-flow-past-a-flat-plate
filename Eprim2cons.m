function [E1, E2, E3, E5] = Eprim2cons(rho, u, p, v, T, mu, lambda, k, c_v, dx, dy, call_case)

% Preallocate other arrays with zeros 
[numy,numx] = size(u);
E1 = zeros(numy,numx);
E2 = zeros(numy,numx);
E3 = zeros(numy,numx);
E5 = zeros(numy,numx);


E1 = rho.*u; % From continuity

tau_xx = xxshear(u, v, lambda, mu, dx, dy, call_case);
E2 = rho.*u.^2 + p - tau_xx; % From x-momentum

tau_xy = xyshear(u, v, mu, dx, dy, call_case);
E3 = rho.*u.*v - tau_xy; % From y-momentum

q_x = heat_flux_x(T, k, dx, call_case);
E5 = (rho.*(c_v*T + (u.^2 + v.^2)/2) + p).*u - u.*tau_xx - v.*tau_xy + q_x; % From energy