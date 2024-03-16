function tau_xy = xyshear(u, v, mu, dx, dy, call_case)
% The xy component of the shear stress

% Preallocate other arrays with zeros 
[numy, numx] = size(u);
du_dy = zeros(numy, numx);
dv_dx = zeros(numy, numx);

if (strcmp(call_case, 'E_predictor') || strcmp(call_case, 'E_corrector'))

    du_dy = ddy_central_updated(u,dy);
    
    if (strcmp(call_case, 'E_predictor'))
        
    dv_dx = ddx_bwd_updated(v,dx);%ddx_bwd_updated(v,dx);
    else
    dv_dx = ddx_fwd_updated(v,dx);%ddx_fwd_updated(v,dx);

    end
    
elseif (strcmp(call_case, 'F_predictor') || strcmp(call_case, 'F_corrector'))
    
     dv_dx = ddx_central_updated(v,dx);


    if (strcmp(call_case, 'F_predictor'))
        
     du_dy = ddy_bwd_updated(u,dy);

    else
        
     du_dy = ddy_fwd_updated(u,dy);

    end
end

tau_xy = mu.*(du_dy + dv_dx);