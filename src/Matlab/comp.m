function u0_comp = comp(u0_proj,Nx,Ny,Lx,Ly)

dx = Lx/Nx; dy = Ly/Ny;

u0_comp = zeros(Nx,Ny);
for iy = 1:Ny
    yj = dy*(iy-1/2);
    for ix = 1:Nx
        xi = dx*(ix-1/2);
        
        for i = 0:Nx-1
            hi = cos(pi*i/Lx*xi);
            for j = 0:Ny-1
                hj = cos(pi*j/Ly*yj);
                u0_comp(ix,iy) = u0_comp(ix,iy) + u0_proj(i+1,j+1)*hi*hj;
            end
        end
    end
end