function int = integral_ij(u0grid,i,j,Nx,Ny,Lx,Ly)

dx = Lx/Nx; dy = Ly/Ny;

int = 0;
for iy = 1:Ny
    yj = dy*(iy-1/2);
    hj = sqrt(1+sign(j))*cos(pi*j/Ly*yj);
                
    for ix = 1:Nx
        xi = dx*(ix-1/2);

        hi = sqrt(1+sign(i))*cos(pi*i/Lx*xi);

        integrand = u0grid(ix,iy)*hi*hj;

        int = int + dx*dy*integrand;
        
    end
end
