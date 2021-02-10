%Corrector step for velocity calculation involving pressure
function [u,v] = corrector(u,v,utemp,vtemp,dt,rho,p,hx,hy,nx,ny,F)

for j = 2:ny+1
    for i = 3:nx+1
        u(i,j) = utemp(i,j) - dt/rho *(p(i,j)-p(i-1,j))*hx +F*dt;
    end
end
for j = 3:ny+1
    for i = 2:nx+1
        v(i,j) = vtemp(i,j) - dt/rho *(p(i,j)-p(i,j-1))*hy;
    end
end

% WALL BC
u(ny+1,:)=0;
u(1,:)=0;
v(:,1)=0;
v(:,ny+1)=0;

%boundary BC
u(:,nx) = u(:,2);
u(:,1) = u(:,nx-1);
v(:,nx) = v(:,2);
v(:,1) = v(:,nx-1);

u(:,1)   =u(:,2)   -2*(u(:,2));
u(:,ny+2)=u(:,ny+1)-2*(u(:,ny+1));
v(1,:)   =v(2,:)   -2*(v(2,:));
v(nx+2,:)=v(nx+1,:)-2*(v(nx+1,:));        


end

