% Create Laplacian operator for solving pressure Poisson equation

function [L] = Laplacian(L,nx,ny,hx,hy)

for j=2:ny
    for i=2:nx
        L(i+(j-1)*nx, i+(j-1)*nx)=2*hx^2+2*hy^2 ;
        for ii=i-1:2:i+1
            if (ii>0 && ii<=nx) %Interior point
                L(i+(j-1)*nx, ii+(j-1)*nx)=-hx^2;
            else %Neuman conditions on boundary
                L(i+(j-1)*nx, i+(j-1)*nx)= L(i+(j-1)*nx, i+(j-1)*nx)-hx^2;
            end
        end
        for jj=j-1:2:j+1
            if (jj>0 && jj<=ny)  %Interior point
                L(i+(j-1)*nx, i+(jj-1)*nx)=-hy^2;
            else %Neuman conditions on boundary
                L(i+(j-1)*nx, i+(j-1)*nx)= L(i+(j-1)*nx, i+(j-1)*nx)-hy^2;
            end
        end
    end
end

end

