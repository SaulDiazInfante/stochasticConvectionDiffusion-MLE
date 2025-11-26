function int = integral_ijkl(v1grid, v2grid, i, j, k, l, Nx, Ny, Lx, Ly)
dx = Lx/Nx; dy = Ly/Ny;

int = 0;
for iy = 1:Ny
    yj = dy * (iy - 1/2);
    hj = cos(pi * j / Ly * yj);
    hj_prime = -pi * j / Ly * sin(pi*j / Ly * yj);
    hl = cos(pi * l / Ly * yj);
                
    for ix = 1:Nx
        xi = dx * (ix - 1 / 2);
        hi = cos(pi * i / Lx * xi);
        hi_prime = -pi*i/Lx*sin(pi*i/Lx*xi);
        hk = cos(pi*k/Lx*xi);

        integrand = v1grid(ix,iy)*hi_prime*hj*hk*hl+v2grid(ix,iy)*hi*hj_prime*hk*hl;
        int = int + dx*dy*integrand;
        
    end
end