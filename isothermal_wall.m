

function [rho, u, v, p, T] = isothermal_wall(rho, u, v, p, T, rho_infinity, u_infinity, p_infinity, T_infinity, R)
% Apply boundary conditions (constant-temperature wall case).

[numy, numx] = size(rho);

% Leading edge 
u(1,1) = 0;
v(1,1) = 0;
T(1,1) = T_infinity;
p(1,1) = p_infinity;
rho(1,1) = rho_infinity;
% 
% % Inflow % comment all inflow
% u(2:numy,1) = u_infinity;
% p(2:numy,1) = p_infinity;
% T(2:numy,1) = T_infinity;
% rho(2:numy,1) = rho_infinity;


% For the upper boundary:
u(numy,2:numx) = u_infinity;
p(numy,2:numx) = p_infinity;
T(numy,2:numx) = T_infinity;
rho(numy,2:numx) = rho_infinity;


% % Outflow %comment all otflow
% u(2:numy-1,numx) = 2*u(2:numy-1,numx-1) - u(2:numy-1,numx-2);
% v(2:numy-1,numx) = 2*v(2:numy-1,numx-1) - v(2:numy-1,numx-2);
% p(2:numy-1,numx) = 2*p(2:numy-1,numx-1) - p(2:numy-1,numx-2);
% T(2:numy-1,numx) = 2*T(2:numy-1,numx-1) - T(2:numy-1,numx-2);
% rho(2:numy-1,numx) = 2*rho(2:numy-1,numx-1) - rho(2:numy-1,numx-2);


% Wall 
 u(1,2:numx) = 0;
 v(1,2:numx) = 0;
 T(1,2:numx) = T_infinity;
 p(1,2:numx) = 2*p(2,2:numx) - p(3,2:numx);
 rho(1,2:numx) = p(1,2:numx)./(R*T(1,2:numx));