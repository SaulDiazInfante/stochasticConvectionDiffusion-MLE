function int = integral_ijkl(v1grid,v2grid,i,j,k,l,Nx,Ny,Lx,Ly)

dx = Lx/Nx; dy = Ly/Ny;

int = 0;
for iy = 1:Ny
    yj = dy*(iy-1/2);
    hj = sqrt(1+sign(j))*cos(pi*j/Ly*yj);        
    %hj_prime = -pi*j/Ly*sqrt(1+sign(j))*sin(pi*j/Ly*yj); %Corrected
    hl = sqrt(1+sign(l))*cos(pi*l/Ly*yj);
    hl_prime = -pi*l/Ly*sqrt(1+sign(l))*sin(pi*l/Ly*yj);

    for ix = 1:Nx
        xi = dx*(ix-1/2);

        hi = sqrt(1+sign(i))*cos(pi*i/Lx*xi);
        %hi_prime = -pi*i/Lx*sqrt(1+sign(i))*sin(pi*i/Lx*xi); %Corrected
        hk = sqrt(1+sign(k))*cos(pi*k/Lx*xi);
        hk_prime = -pi*k/Lx*sqrt(1+sign(k))*sin(pi*k/Lx*xi);

        %integrand = v1grid(ix,iy)*hi_prime*hj*hk*hl+v2grid(ix,iy)*hi*hj_prime*hk*hl;
        integrand = v1grid(ix,iy)*hi*hj*hk_prime*hl+v2grid(ix,iy)*hi*hj*hk*hl_prime;
        int = int + dx*dy*integrand;
        
    end
end