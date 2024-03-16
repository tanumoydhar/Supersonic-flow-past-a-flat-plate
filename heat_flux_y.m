

function q_y = heat_flux_y(T, k, dy, call_case)
% Compute the y-component of the heat flux

% Preallocate other arrays with zeros 
[numy, numx] = size(T);
dT_dy = zeros(numy, numx);

if (strcmp(call_case, 'F_predictor'))
    
dT_dy = ddy_bwd_updated(T,dy);

elseif (strcmp(call_case, 'F_corrector'))
    
dT_dy = ddy_fwd_updated(T,dy);

end

q_y = - k.*dT_dy;