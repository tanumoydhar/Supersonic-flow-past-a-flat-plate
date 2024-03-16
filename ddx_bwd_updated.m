function dfdx = ddx_bwd_updated(f,dx)
    
    % determine field size
    [ny,nx]     = size(f);

    % allocate return field
    dfdx        = zeros(ny,nx);
    
    % backward difference
    for i=2:nx
        for j=1:ny
            dfdx(j,i) = (f(j,i)-f(j,i-1))/dx;
        end
    end

%     % forward difference for first point
%     i = 1;
%     for j=1:ny
%         dfdx(j,i) = (f(j,i+1)-f(j,i))/dx;
%     end
    
    %     % assuming periodicity  (left boudary)
        i = 1;
        for j=1:ny
            dfdx(j,i) = (f(j,i)-f(j,end))/dx;
        end
end