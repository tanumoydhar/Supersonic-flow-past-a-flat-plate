clear all; close all; clc;

%Tanumoy Dhar, MAE 290C 
% 2-D compressible Navier-Stokes equations with closure for calorically perfect gas using MacCormack's method. 
%  the problem of the flow over a flat plate with a sharp leading-edge at Mach 4. 


%Question

%Your code must have the following features:

% During runtime, update a plot that shows the fields of ρ,u,v,e,p and T in subplots. 
% Rendering graphical output will slow down your code (contourf in particular). You might 
% want to consider updating the figure only every 10 or 50 or so time steps (mod is your best
% friend here). Also, add a line plot that tracks the convergence of the flow field towards
% the steady-state solution. You may pick a quantity of your choice. Be
% sure to use labels, colorbars, etc. - Done!


% Consistent 2nd-order of accuracy discretization everywhere in the flow
% field and at the boundaries  - Done!


% Clear, readable, and well documented. Use lots of comments - Done!


% Make smart use of self-written functions that you may store in external m-files.
% In particular, you are allowed and encouraged to use your self-written functions
% from previous assignments! 


% Efficiency is not the highest priority—we prefer a clean code. 
% That does not mean that you can be wasteful. Don't have your code do any unnecessary operations 
% that severely affect its time-to-solution. For example, computing values at the boundaries that 
% your code doesn't really need (because they get overwritten by boundary conditions) is OK if it
% helps clarity, but computing entire flow fields that are not used is
% wasteful  - Done!



% Define physical parameters and simulation parameters

numx = 75;                               % Number of grid points along the x-axis
numy = 80;                               % Number of grid points along the y-axis
length_flat_plate = 1e-5;                % Flat plate length [m]
height_flat_plate = 8e-6;                % Height for the computational domain [m]
dx = length_flat_plate/(numx - 1);       %  grid points spacing along the x-coordinate [m]
dy = height_flat_plate/(numy - 1);       %  grid points spacing along the y-coordinate [m]
x = 0:dx:length_flat_plate;              % x-coordinate of the domain
y = 0:dy:height_flat_plate;              % y-coordinate of the domain
total_time_steps = 1500;                 % number of time steps  
[xx,yy] = ndgrid(x,y);
xlim = ([min(min(x)), max(max(x))]);
ylim = ([min(min(y)), max(max(y))]);
t = 1;    % Time step counter
time = 0; % Physical time


gamma = 1.4;  mu_0 = 1.735e-5; 
T_0 = 288.15;  Pr = 0.71; R = 287; 

% Compute initial physical parameters (speed of sound, viscosity, etc)

Ma_freestream = 4.00;                   % Freestream Mach number
sndspd_freestream = sqrt(gamma*R*T_0); % Freestream speed of sound [m/s]
pressure_freestream = 101300;          % Freestream pressure [Pa]
temp_freestream = T_0;                 % Freestream temperature [K]


c_v = R/(gamma - 1);  c_p = gamma*c_v;    
rho_inf = pressure_freestream/(R*temp_freestream); 
mu = mu_0*(temp_freestream/T_0).^(3/2).*((T_0 + 110.4)./(temp_freestream + 110.4));
Re_L = rho_inf*Ma_freestream*sndspd_freestream*length_flat_plate./mu;


%Initialize grid
p = ones(numy,numx)*pressure_freestream;
rho = ones(numy,numx)*rho_inf;
T = ones(numy,numx)*temp_freestream;
T(1,:) = temp_freestream;
u = 0.01*((rand(numy,numx)-0.5)*2);%*ones(numy,numx)

u(1,:) = 0;                              % No-slip boundary condition 

v = zeros(numy,numx);
mu = mu_0*(T/T_0).^(3/2)*(T_0 + 110.4)./(T + 110.4);
lambda = - 2/3*mu;                       % Second viscosity
k = mu*(c_p/Pr);                           % thermal conductivity


%Initialize  variables

U1_p  = zeros(numy, numx);
U2_p  = zeros(numy, numx);
U3_p  = zeros(numy, numx);
U5_p  = zeros(numy, numx);
rho_p = zeros(numy, numx);
u_p   = zeros(numy, numx);
v_p   = zeros(numy, numx);
T_p   = zeros(numy, numx);    
p_p   = zeros(numy, numx);


while t < total_time_steps
    % Make sure conservative and primitive variables are up to date

    delta_t = 2.35e-11;
    time = time + delta_t;
    t = t + 1;

    % Predictor %*********************************************************
    % Update all necessary physical parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     [U1, U2, U3, U5] = Uprim2cons(rho, u, v, T, c_v);
     [E1, E2, E3, E5] = Eprim2cons(rho, u, p, v, T, mu, lambda, k, c_v, dx, dy, 'E_predictor');
     [F1, F2, F3, F5] = Fprim2cons(rho, u, p, v, T, mu, lambda, k, c_v, dx, dy, 'F_predictor');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     U1_p = U1 - delta_t*(ddx_fwd_updated(E1,dx) + ddy_fwd_updated(F1,dy));
     U2_p = U2 - delta_t*(ddx_fwd_updated(E2,dx) + ddy_fwd_updated(F2,dy));
     U3_p = U3 - delta_t*(ddx_fwd_updated(E3,dx) + ddy_fwd_updated(F3,dy));
     U5_p = U5 - delta_t*(ddx_fwd_updated(E5,dx) + ddy_fwd_updated(F5,dy));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
    % Update primitive variables
    
    [rho_p(2:numy-1,2:numx-1), u_p(2:numy-1,2:numx-1), v_p(2:numy-1,2:numx-1),...
    T_p(2:numy-1,2:numx-1)] = Ucons2prim(U1_p(2:numy-1,2:numx-1), ...
    U2_p(2:numy-1,2:numx-1), U3_p(2:numy-1,2:numx-1), U5_p(2:numy-1,2:numx-1), c_v);
    
    p_p(2:numy-1,2:numx-1) = rho_p(2:numy-1,2:numx-1)*R.*T_p(2:numy-1,2:numx-1);
    
    
    % Enforce BCs on primitive variables
    [rho_p, u_p, v_p, p_p, T_p] = isothermal_wall(rho_p, u_p, v_p, p_p, T_p, rho_inf,...
    Ma_freestream*sndspd_freestream, pressure_freestream, temp_freestream, R);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    mu_p = mu_0*(T_p/T_0).^(3/2)*(T_0 + 110.4)./(T_p + 110.4);
    lambda_p = - 2/3*mu_p;
    k_p = mu_p*c_p/Pr;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Corrector %%******************************************************
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        
   [E1_p, E2_p, E3_p, E5_p] = Eprim2cons(rho_p, u_p, p_p, v_p, T_p, mu_p, ...
       lambda_p, k_p, c_v, dx, dy, 'E_corrector');
   [F1_p, F2_p, F3_p, F5_p] = Fprim2cons(rho_p, u_p, p_p, v_p, T_p, mu_p, ...
       lambda_p, k_p, c_v, dx, dy, 'F_corrector');

       
     U1 = 0.5*(U1_p+U1) - delta_t*(ddx_bwd_updated(E1_p,dx) + ddy_bwd_updated(F1_p,dy));
     U2 = 0.5*(U2_p+U2) - delta_t*(ddx_bwd_updated(E2_p,dx) + ddy_bwd_updated(F2_p,dy));
     U3 = 0.5*(U3_p+U3) - delta_t*(ddx_bwd_updated(E3_p,dx) + ddy_bwd_updated(F3_p,dy));
     U5 = 0.5*(U5_p+U5) - delta_t*(ddx_bwd_updated(E5_p,dx) + ddy_bwd_updated(F5_p,dy));
    
    % Update primitive variables

    [rho(2:numy-1,2:numx-1), u(2:numy-1,2:numx-1), v(2:numy-1,2:numx-1),...
        T(2:numy-1,2:numx-1)] =...
        Ucons2prim(U1(2:numy-1,2:numx-1), ...
        U2(2:numy-1,2:numx-1), U3(2:numy-1,2:numx-1), U5(2:numy-1,2:numx-1), c_v);
    
    
    p(2:numy-1,2:numx-1) = rho(2:numy-1,2:numx-1)*R.*T(2:numy-1,2:numx-1);
    
    % Enforce BCs on primitive variables

    [rho, u, v, p, T] = isothermal_wall(rho, u, v, p, T, rho_inf,...
        Ma_freestream*sndspd_freestream, pressure_freestream, temp_freestream, R);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    
    mu = mu_0*(T/T_0).^(3/2)*(T_0 + 110.4)./(T + 110.4);
    lambda = - 2/3*mu;
    k = mu*(c_p/Pr);    
       
    
    et = c_v.*T + 0.5*(u.^2 + v.^2);
    local_sp_speed = sqrt(gamma*R*T); 
    M = sqrt(u.^2 + v.^2)./local_sp_speed;
  
 % Data output/visualization
   
    
    
     if mod(t,15) == 0 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,3,1)
    [xx,yy] = ndgrid(x,y);
    GTplot = pcolor(xx,yy,rho');
    set(GTplot, 'EdgeColor', 'none');
    colormap(jet(1024))
    shading interp
    axis equal tight
    colorbar()
    xlabel('$x$','fontsize',15,'interpreter','latex');
    ylabel('$y$','fontsize',15,'interpreter','latex');
    title('$\rho[kg/m^3]$','fontsize',10,'interpreter','latex');
     


    subplot(2,3,2)
    [xx,yy] = ndgrid(x,y);
    GTplot = pcolor(xx,yy,u');
    set(GTplot, 'EdgeColor', 'none');
    colormap(jet(1024))
    shading interp
    axis equal tight
    colorbar()
    xlabel('$x$','fontsize',15,'interpreter','latex');
    ylabel('$y$','fontsize',15,'interpreter','latex');
    title('$u[m/s]$','fontsize',10,'interpreter','latex');


    subplot(2,3,3)
    [xx,yy] = ndgrid(x,y);
    GTplot = pcolor(xx,yy,v');
    set(GTplot, 'EdgeColor', 'none');
    colormap(jet(1024))
    shading interp
    axis equal tight
    colorbar()
    xlabel('$x$','fontsize',15,'interpreter','latex');
    ylabel('$y$','fontsize',15,'interpreter','latex');
    title('$v[m/s]$','fontsize',10,'interpreter','latex');


    subplot(2,3,4)
    [xx,yy] = ndgrid(x,y);
    GTplot = pcolor(xx,yy,et');
    set(GTplot, 'EdgeColor', 'none');
    colormap(jet(1024))
    shading interp
    axis equal tight
    colorbar()
    xlabel('$x$','fontsize',15,'interpreter','latex');
    ylabel('$y$','fontsize',15,'interpreter','latex');
    title('$e[J/kg]$','fontsize',10,'interpreter','latex');


    subplot(2,3,5)
    [xx,yy] = ndgrid(x,y);
    GTplot = pcolor(xx,yy,p');
    set(GTplot, 'EdgeColor', 'none');
    colormap(jet(1024))
    shading interp
    axis equal tight
    colorbar()
    xlabel('$x$','fontsize',15,'interpreter','latex');
    ylabel('$y$','fontsize',15,'interpreter','latex');
    title('$P[Pa]$','fontsize',10,'interpreter','latex');


    subplot(2,3,6)
    [xx,yy] = ndgrid(x,y);
    GTplot1 = pcolor(xx,yy,T');
    set(GTplot1, 'EdgeColor', 'none');
    colormap(jet(1024))
    shading interp
    axis equal tight
    colorbar()
    xlabel('$x$','fontsize',15,'interpreter','latex');
    ylabel('$y$','fontsize',15,'interpreter','latex');
    title('$T[K]$','fontsize',10,'interpreter','latex');
    
    drawnow
     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end



%% Solution at 1500 time steps
subplot(2,3,1)
[xx,yy] = ndgrid(x,y);
GTplot = pcolor(real(xx),real(yy),real(rho)');
set(GTplot, 'EdgeColor', 'none');
colormap(jet(1024))
shading interp
axis equal tight
colorbar()
xlabel('$x$','fontsize',15,'interpreter','latex');
ylabel('$y$','fontsize',15,'interpreter','latex');
title('$\rho[kg/m^3]$','fontsize',10,'interpreter','latex');


subplot(2,3,2)
[xx,yy] = ndgrid(x,y);
GTplot = pcolor(xx,yy,u');
set(GTplot, 'EdgeColor', 'none');
colormap(jet(1024))
shading interp
axis equal tight
colorbar()
xlabel('$x$','fontsize',15,'interpreter','latex');
ylabel('$y$','fontsize',15,'interpreter','latex');
title('$u[m/s]$','fontsize',10,'interpreter','latex');


subplot(2,3,3)
[xx,yy] = ndgrid(x,y);
GTplot = pcolor(xx,yy,v');
set(GTplot, 'EdgeColor', 'none');
colormap(jet(1024))
shading interp
axis equal tight
colorbar()
xlabel('$x$','fontsize',15,'interpreter','latex');
ylabel('$y$','fontsize',15,'interpreter','latex');
title('$v[m/s]$','fontsize',10,'interpreter','latex');


subplot(2,3,4)
[xx,yy] = ndgrid(x,y);
GTplot = pcolor(xx,yy,et');
set(GTplot, 'EdgeColor', 'none');
colormap(jet(1024))
shading interp
axis equal tight
colorbar()
xlabel('$x$','fontsize',15,'interpreter','latex');
ylabel('$y$','fontsize',15,'interpreter','latex');
title('$e[J/kg]$','fontsize',10,'interpreter','latex');


subplot(2,3,5)
[xx,yy] = ndgrid(x,y);
GTplot = pcolor(xx,yy,p');
set(GTplot, 'EdgeColor', 'none');
colormap(jet(1024))
shading interp
axis equal tight
colorbar()
xlabel('$x$','fontsize',15,'interpreter','latex');
ylabel('$y$','fontsize',15,'interpreter','latex');
title('$P[Pa]$','fontsize',10,'interpreter','latex');


subplot(2,3,6)
[xx,yy] = ndgrid(x,y);
GTplot1 = pcolor(xx,yy,T');
set(GTplot1, 'EdgeColor', 'none');
colormap(jet(1024))
shading interp
axis equal tight
colorbar()
xlabel('$x$','fontsize',15,'interpreter','latex');
ylabel('$y$','fontsize',15,'interpreter','latex');
title('$T[K]$','fontsize',10,'interpreter','latex');


%% Numerical Schlieren 

%Schlieren photography leverages variations in the refractive index caused
% by density gradients. It is a widely used experimental technique to visualize 
% shocks. Numerical schlieren refers to fake schlieren images generated from 
% numerical data. 


figure()
grad_rho = sqrt((ddx_fwd_updated(rho,dx)).^2 + (ddy_fwd_updated(rho,dy)).^2);
beta = 0.8; kappa = 10;
sch = beta.*exp(-kappa.*(grad_rho./max(max(grad_rho))));
GTplot = pcolor(xx,yy,sch');
set(GTplot, 'EdgeColor', 'none');
colormap(gray(1024))
shading interp
axis equal tight
colorbar()
xlabel('$x$','fontsize',25,'interpreter','latex');
ylabel('$y$','fontsize',25,'interpreter','latex');
hold on

theta = asin(1/Ma_freestream); %using the Ma_freestream
x1 = linspace(0,length_flat_plate,1000);
y1 = tan(theta).*x1;
plot(x1,y1,'-.', 'linewidth', 1.50)





%%
% Implement the adiabatic wall boundary condition into your code. Make sure
% that 2nd-order accuracy is maintained. Run the Mach 4 case for an adiabatic wall and

% Compare the normalized temperature and pressure profiles (as a function of
% y) for the two boundary conditions at three different x-locations,xL=0.25,0.50,0.75. 
% Compare the wall temperature (as a function of x) for both cases in a second plot.








