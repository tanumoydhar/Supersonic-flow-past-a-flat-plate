function tau_xx = xxshear(u, v, lambda, mu, dx, dy, call_case)
% Calculate the x-direction of the normal stress.

% Preallocate other arrays with zeros 
[numy, numx] = size(u);
du_dx = zeros(numy, numx);
dv_dy = zeros(numy, numx);

if (strcmp(call_case, 'E_predictor'))
    
du_dx = ddx_bwd_updated(u,dx);%ddx_bwd_updated(u,dx);

elseif (strcmp(call_case, 'E_corrector'))
    
du_dx = ddx_fwd_updated(u,dx);%ddx_fwd_updated(u,dx); 

end

dv_dy = ddy_central_updated(v,dy);

tau_xx = lambda.*(du_dx + dv_dy) + 2*mu.*du_dx;