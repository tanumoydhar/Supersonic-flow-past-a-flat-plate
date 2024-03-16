function dfdy = ddy_bwd_updated(f,dy)
 
    % determine field size
    [ny,nx]     = size(f);

    % allocate return field
    dfdy        = zeros(ny,nx);
    
    % forward difference
    for i=1:nx
        for j=2:ny
            dfdy(j,i) = (f(j,i)-f(j-1,i))/dy;
        end
    end
    
    % backward difference for last point
    j = 1;
    for i=1:nx
        dfdy(j,i) = (f(j+1,i)-f(j,i))/dy;
    end
    
%     % assuming periodicity (right boudary)
%     i = nx;
%     for j=1:ny
%         dfdx(i,j) = (f(1,j)-f(i,j))/dx;
%     end
end




