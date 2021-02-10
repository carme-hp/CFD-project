%This code solves the incompressible unsteady 2D Navier-Stokes equation on
% a uniform staggered grid --> pressure is stored at the cell center and
% velocities on the faces; this allows for the solution to have a tight
% coupling
%inlet on left side, outlet on right side, top and bottom are walls

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

% Note that when viewed separately, the p-grid is a NxN while the u-grid
% is (N-1)x N and the v-grid is N x(N-1)
% dx is the spacing of each grid
% the entire domain size will be N+2, the inner nodes are N (this is due
% to the staggered grid setup, we need an additional boundary node for the
% velocities)

%STEPS________________________________
% Set input parameters: viscosity, density, number of grid points, 
% time information, and boundary conditions
% Create the index extents and the computational grid 
% Initialize any arrays you use to allocate the memory
% Create the Laplacian operator 
% Loop over time 
% Apply boundary conditions to the velocity 
% Perform the predictor step 
% Form the right-hand-side of the Poisson equation 
% Solve for the pressure using \
% Perform the corrector step 
% Plot  
% End Simulation

clear; close all; clc;

%% set input parameters
nu = 0.001; %kinematic viscosity
rho = 1;  %density
nx = 20;  % node points excluding boundaries
ny = 20;
Lx = 2;   %domain length over x
Ly = 2;   %domain height over y
dt = 0.01;   % timestep
time = 10; %end time (s)
actualtime = 0;
F=10;
%% define domain & generate grid

x (2:nx+2)=linspace(0,Lx,nx+1);
y (2:ny+2)=linspace(0,Ly,ny+1);
% Define location of the middle points for pressure values
xm(2:nx+1)=0.5*(x(2:nx+1) + x(3:nx+2));
ym(2:ny+1)=0.5*(y(2:ny+1) + y(3:ny+2));
% Create mesh sizes
dx=x(3)-x(2);
dy=y(3)-y(2);
hx=1/dx;
hy=1/dy;
%% Memory allocation to matrices
utemp = zeros(nx+2,ny+2);
vtemp = zeros(nx+2,ny+2);
u = zeros(nx+2,ny+2);
v = zeros(nx+2,ny+2);
p = zeros(nx+1,ny+1);
R = zeros(nx+1,1);
L=zeros(nx*ny,nx*ny);

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
        % temporary velocity calculation (predictor step)
        [utemp,vtemp]=predictor(utemp,vtemp,u,v,nx,ny,hx,hy,nu,dt,F);
        
        
        %% compute RHS
        n = 0;
        for j = 2:ny+1
            for i = 2:nx+1
                n = n+1;
                R(n) = -rho/dt*((utemp(i+1,j)-utemp(i,j))*hx ...
                    +(vtemp(i,j+1)-vtemp(i,j))*hy);
            end
        end
        %% solve pressure
        % pressure solve with Laplacian pv = L\RHS
        % pv is a vector including the pressure value of every cell
        pv = L\R;
        %convert pv into pressure field
        n = 0;
        for j = 2:ny+1
            for i = 2:nx+1
                n = n+1;
                p(i,j) = pv(n);
            end
        end
        
        p(:,ny+1) =p(:,ny);	
        p(:,1) = p(:,2);
        p(:,nx) = p(:,2);
        p(:,1) = p(:,nx-1);
        
        %% velocity correction (corrector step)
        [u,v] = corrector(u,v,utemp,vtemp,dt,rho,p,hx,hy,nx,ny,F);
                
        %% plot
        figure(1); clf(1);
        title('Velocity & Pressure Field','fontsize',20);
        xlim([0,Lx]);
        ylim([0,Ly]);
        hold on
        %pressure plot
        %contourf(x(2:nx+1),y(2:ny+1),p(2:nx+1,2:ny+1)');
        %velocity plot
        contourf(x(2:nx+1),y(2:ny+1),u(2:nx+1,2:ny+1)');
        %velocity vectors
        quiver(x(2:nx+1),y(2:ny+1),...
            u(2:nx+1,2:ny+1)',v(2:nx+1,2:ny+1)','filled', 'k');
        xmax=Lx-dx;
        ymax=Ly-dy;
        % Set Limits on X & Y Axis
        xlim([0 xmax]);
        ylim([0 ymax]);
        col = colorbar;
        %col.Label.String = 'Pressure (Pa)';
        drawnow
        pause(0.01)
        %% timestep end
        actualtime = actualtime + dt;
    end
end % delta t criterion check end
