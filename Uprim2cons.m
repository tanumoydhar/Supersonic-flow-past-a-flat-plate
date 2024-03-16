
function [U1, U2, U3, U5] = Uprim2cons(rho, u, v, T, c_v)
% Calculate the flux vector F from the primitive flow field variables

% Preallocate other arrays with zeros 
[numy,numx] = size(u);
U1 = zeros(numy,numx);
U2 = zeros(numy,numx);
U3 = zeros(numy,numx);
U5 = zeros(numy,numx);

U1 = rho; % From continuity
U2 = rho.*u; % From momentum-x
U3 = rho.*v; % From momentum-y
U5 = rho.*(c_v*T + (u.^2 + v.^2)/2); % From energy