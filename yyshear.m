

function tau_yy = yyshear(u, v, lambda, mu, dx, dy, call_case)
% Compute the y-component of the normal stress.

% Preallocate other arrays with zeros 
[numy, numx] = size(v);
du_dx = zeros(numy, numx);
dv_dy = zeros(numy, numx);


if (strcmp(call_case, 'F_predictor'))
dv_dy = ddy_bwd_updated(v,dy);


elseif (strcmp(call_case, 'F_corrector'))
dv_dy = ddy_fwd_updated(v,dy);

end

du_dx = ddx_central_updated(v,dx);

tau_yy = lambda.*(du_dx + dv_dy) + 2.0.*mu.*dv_dy;