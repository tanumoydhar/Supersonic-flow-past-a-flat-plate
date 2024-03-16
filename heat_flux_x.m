function q_x = heat_flux_x(T, k, dx, call_case)
% Compute the x-component of the heat flux

% Preallocate other arrays with zeros 
[numy, numx] = size(T);
dT_dx = zeros(numy, numx);

if (strcmp(call_case, 'E_predictor'))
    
 dT_dx = ddx_bwd_updated(T,dx);
 
elseif (strcmp(call_case, 'E_corrector'))
    
 dT_dx = ddx_fwd_updated(T,dx);
 
end

q_x = - k.*dT_dx;