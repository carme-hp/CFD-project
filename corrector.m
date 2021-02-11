%Corrector step for velocity calculation involving pressure
function [u,v] = corrector(u,v,utemp,vtemp,dt,rho,p,hx,hy,nx,ny)

for j = 2:ny
    for i = 2:nx
        u(i,j) = utemp(i,j) - dt/rho *(p(i,j)-p(i-1,j))*hx;
    end
end
for j = 2:ny
    for i = 2:nx
        v(i,j) = vtemp(i,j) - dt/rho *(p(i,j)-p(i,j-1))*hy;
    end
end

end

