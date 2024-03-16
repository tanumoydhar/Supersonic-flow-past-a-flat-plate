function dfdy = ddy_central_updated(f,dy)

    % determine field size
    [ny,nx]     = size(f);

    % allocate return field
    dfdy        = zeros(ny,nx);
    
    % central difference
    for i=1:nx
        for j=2:ny-1
            dfdy(j,i) = (f(j+1,i)-f(j-1,i))/2/dy;
        end
    end
    
    % forward difference for first point
    j = 1;
    for i=1:nx
        dfdy(j,i) = (-3*f(j,i)+4*f(j+1,i)-f(j+2,i))/2/dy;
    end
    
    % backward difference for last point
    j = ny;
    for i=1:nx
        dfdy(j,i) = (3*f(j,i)-4*f(j-1,i)+f(j-2,i))/2/dy;
    end
    
%     % assuming periodicity (left boudary)
%     i = 1;
%     for j=1:ny
%         dfdx(i,j) = (f(i+1,j)-f(end,j))/2/dx;
%     end
%     % assuming periodicity (right boudary)
%     i = nx;
%     for j=1:ny
%         dfdx(i,j) = (f(1,j)-f(i-1,j))/dx;
%     end
end