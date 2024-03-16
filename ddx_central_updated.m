function dfdx = ddx_central_updated(f,dx)

    % determine field size
    [ny,nx]     = size(f);

    % allocate return field
    dfdx        = zeros(ny,nx);
    
    % central difference
    for i=2:nx-1
        for j=1:ny
            dfdx(j,i) = (f(j,i+1)-f(j,i-1))/2/dx;
        end
    end
    
%     % forward difference for first point
%     i = 1;
%     for j=1:ny
%         dfdx(j,i) = (-3*f(j,i)+4*f(j,i+1)-f(j,i+2))/2/dx;
%     end
%     
%     % backward difference for last point
%     i = nx;
%     for j=1:ny
%         dfdx(j,i) = (3*f(j,i)-4*f(j,i-1)+f(j,i-2))/2/dx;
%     end
    
    % assuming periodicity (left boudary)
    i = 1;
    for j=1:ny
        dfdx(j,i) = (f(j,i+1)-f(j,end))/2/dx;
    end
    % assuming periodicity (right boudary)
    i = nx;
    for j=1:ny
        dfdx(j,i) = (f(j,1)-f(j,i-1))/dx;
    end
end