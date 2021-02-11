%This function calculates velocity prediction on staggered grid

function [utemp,vtemp] = predictor(utemp,vtemp,u,v,nx,ny,hx,hy,nu,dt)
%temporary u %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 2:ny
            for i = 2:nx
                %calculate v values for u's staggered grid
                vstagg = 0.25*(v(i-1,j)+v(i-1,j+1)+v(i,j)+v(i,j+1));
                %substitute into u momentum
                utemp(i,j) = u(i,j) + dt* ...
                    (nu*(u(i-1,j)-2*u(i,j)+u(i+1,j))*hx^2 ...
                    +nu*(u(i,j-1)-2*u(i,j)+u(i,j+1))*hy^2 ...
                    -u(i,j)*(u(i+1,j)-u(i-1,j))*0.5*hx ...
                    -vstagg*(u(i,j+1)-u(i,j-1))*0.5*hy);
            end
        end
%temporary v %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 2:ny
            for i = 2:nx
                %calculate u values for v's staggered grid
                ustagg = 0.25*(u(i,j-1)+u(i,j)+u(i+1,j-1)+u(i+1,j));
                %substitute into v momentum
                vtemp(i,j) = v(i,j)+dt* ...
                    (nu*(v(i-1,j)-2*v(i,j)+v(i+1,j))*hx^2 ...
                    +nu*(v(i,j-1)-2*v(i,j)+v(i,j+1))*hy^2 ...
                    -ustagg*(v(i+1,j)-v(i-1,j))*0.5*hx ...
                    -v(i,j)*(v(i,j+1)-v(i,j-1))*0.5*hy);
            end
        end
end

