% This code solves the incompressible transient 2D Navier-Stokes equation on
% a uniform staggered grid --> pressure is stored at the cell center and
% velocities on the faces; this allows for the solution to have a tight
% coupling

%   staggered grid
%--------------------------------------
%   -p- -u- -p- -u- -p- -u- -p- -u- -p-
%   -v-     -v-     -v-     -v-     -v-
%   -p- -u- -p- -u- -p- -u- -p- -u- -p-
%   -v-     -v-     -v-     -v-     -v-
%   -p- -u- -p- -u- -p- -u- -p- -u- -p-
%   -v-     -v-     -v-     -v-     -v-
%   -p- -u- -p- -u- -p- -u- -p- -u- -p-
%   -v-     -v-     -v-     -v-     -v-
%   -p- -u- -p- -u- -p- -u- -p- -u- -p-
%---------------------------------------
clear; close all; clc;

%% set input parameters
nu = 0.001; %kinematic viscosity
rho = 1;  %density
nx = 20;  % node points excluding boundaries
ny = 20;
Lx = 2;   %domain length over x
Ly = 2;   %domain height over y
dt = 0.01;   % timestep
time = 1; %end time (s)
actualtime = 0; 
% Create mesh sizes
hx=Lx/nx;
hy=Ly/ny;
%% Memory allocation to matrices
utemp = zeros(nx+1,ny+1); % CARME: nx+2, ny+2 (thats what you get the index error (holds for u.v,p)
vtemp = zeros(nx+1,ny+1);
u = zeros(nx+1,ny+1);
v = zeros(nx+1,ny+1);
p = zeros(nx+1,ny+1);
R = zeros(nx,1); % CARME: wrong dimensions (nx*ny)
L=zeros(nx*ny,nx*ny); % CARME : unnecessary line

%% Laplacian operator
L = Laplacian(L,nx,ny,hx,hy);
% Set pressure BC
L(1,:)=0; L(1,1)=1; 
%% time loop warning for too large timestep
if (max(max(u)))*dt/hx >= 1 || (max(max(v)))*dt/hy >= 1
    disp ('timestep too large, divergence may occur')
else
    %% time loop
    while actualtime<=time
        %% setting of Boundary Conditions
        u(:,1) = 1; %First column is inlet --> not calculated later on
        v(:,1) = 0;
        % BC wall
        u(1,:)=0;
        v(1,:)=0;
        u(nx+1,:)=0; % CARME nx+2
        v(nx+1,:)=0;
        %% temporary velocity calculation (predictor step)
        [utemp,vtemp]=predictor(utemp,vtemp,u,v,nx,ny,hx,hy,nu,dt);
        %Gaussian Layer
        u(:,nx+1) = u(:,nx); % CARME u(:,nx+2) = u(:,nx+1);
        v(:,nx+1) = v(:,nx);
        %% compute RHS
        n = 0;
        for j = 2:ny+1
            for i = 2:nx+1
                n = n+1;
                 R(n) = -rho/dt*((utemp(i+1,j)-utemp(i,j))*hx ...
                     +(vtemp(i,j+1)-vtemp(i,j))*hy);
            end
        end
        R(nx+1,1)=R(nx,1); % CARME:  rhs doesnt include ghost layers (nx*ny) so you dont need that (plus if you did it would be nx+2)
        %% solve pressure
        % pressure solve with Laplacian pv = L\RHS
        % pv is a vector including the pressure value of every cell
        pv = L\R;
        %convert pv into pressure field
        n = 0; % CARME: achtung!! term R(1) corresponds to entry (2,2) in the p-matrix because of ghost layers. 
        for j = 2:ny
            for i = 2:nx
                n = n+1;
                p(i,j) = pv(n);
            end
        end
        p(:,nx+1) = p(:,nx); 
        %% velocity correction (corrector step)
        [u,v] = corrector(u,v,utemp,vtemp,dt,rho,p,hx,hy,nx,ny);
        %Gaussian Layer
        u(:,nx+1) = u(:,nx);
        v(:,nx+1) = v(:,nx);    
        %% plot
        figure(1); clf(1);
        title('Velocity & Pressure Field','fontsize',20);
        xlim([0,Lx]);
        ylim([0,Ly]);
        hold on
        
        x (2:nx+2)=linspace(0,Lx,nx+1);
        y (2:ny+2)=linspace(0,Ly,ny+1);
        
        contourf(x(2:nx+1),y(2:ny+1),p(2:nx+1,2:ny+1)');
        quiver(x(2:nx+1),y(2:ny+1),...
            u(2:nx+1,2:ny+1)',v(2:nx+1,2:ny+1)','filled', 'k');
        xmax=Lx-hx;
        ymax=Ly-hy;
        % Set Limits on X & Y Axis
        xlim([0 xmax]);
        ylim([0 ymax]);
        col = colorbar;
        col.Label.String = 'Pressure (Pa)';
        drawnow
        pause(0.01)
        %% timestep end
        actualtime = actualtime + dt;
    end
end % delta t criterion check end
