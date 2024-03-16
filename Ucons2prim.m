
function [rho, u, v, T] = Ucons2prim(U1, U2, U3, U5, c_v)
% Calculate  the primitive flow field variables from the flux vector U.

rho = U1;   % From continuity
u = U2./U1; % From momentum-x
v = U3./U1; % From momentum-y
T = (U5./U1 - ((U2./U1).^2 + (U3./U1).^2)/2)/c_v;