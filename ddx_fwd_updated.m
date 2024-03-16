function dfdx = ddx_fwd_updated(f,dx)
 
    % determine field size
    [ny,nx]     = size(f);

    % allocate return field
    dfdx        = zeros(ny,nx);
    
    % forward difference
    for i=1:nx-1
        for j=1:ny
            dfdx(j,i) = (f(j,i+1)-f(j,i))/dx;
        end
    end
    
%     % backward difference for last point
%     i = nx;
%     for j=1:ny
%         dfdx(j,i) = (f(j,i)-f(j,i-1))/dx;
%     end
    
%     % assuming periodicity (right boudary)
    i = nx;
    for j=1:ny
        dfdx(j,i) = (f(j,1)-f(j,i))/dx;
    end
end