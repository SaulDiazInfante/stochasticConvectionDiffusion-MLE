% load('Coeficientes.mat')
% vcoef
% stop
mkdir Data

%% Definicion de campo de velocidades:

Lx = 5; Ly = 5;


%% Parametros:

gamma = 2;
Nx = 50;
Ny = 50;

%% Campo de velocidades


xgrid = zeros(Nx,Ny);
ygrid = zeros(Nx,Ny);
v1grid = zeros(Nx,Ny);
v2grid = zeros(Nx,Ny);
u0grid = zeros(Nx,Nx);
for iy = 1:Ny
    for ix = 1:Nx
        xgrid(ix,iy) = Lx / Nx * (ix - 1/2);
        ygrid(ix,iy) = Ly / Ny * (iy - 1/2);
        v1grid(ix,iy) = v1(xgrid(ix, iy), ygrid(ix, iy));
        v2grid(ix,iy) = v2(xgrid(ix, iy), ygrid(ix, iy));
        
        xi = xgrid(ix,iy);
        yj = ygrid(ix,iy);
        u0grid(ix,iy) = u0(xi,yj,Lx,Ly);
    end
end


%% Descomposicion de las condiciones iniciales u0

u0_proj = zeros(Nx,Ny);
for i = 0:Nx-1
    progress_bar(i+1, Nx, 40, 'Decomposing initial condition')
    for j = 0:Ny-1
        u0_proj(i+1,j+1) = (1 / (Lx * Ly)) * integral_ij(u0grid, i, j, Nx, Ny, Lx, Ly);
    end
end

u0_comp_grid = comp(u0_proj,Nx,Ny,Lx,Ly);

figure(1)
subplot(2,2,1)
hold on
title('u_0')
quiver(xgrid,ygrid,v1grid,v2grid)
contour(xgrid,ygrid,u0grid,100); colorbar;
subplot(2,2,2)
hold on
title('u_0 composed')
quiver(xgrid,ygrid,v1grid,v2grid)
contour(xgrid,ygrid,u0_comp_grid,100); colorbar;
axis([0,Lx,0,Ly])
subplot(2,2,3)
hold on
title('Difference')
quiver(xgrid,ygrid,v1grid,v2grid)
contour(xgrid,ygrid,abs(u0_comp_grid-u0grid),100); colorbar;
axis([0,Lx,0,Ly])

clf 

% figure(2)
% contour(u0_proj); colorbar;
% min(min(abs(u0_proj)))
% 
% stop

%% Calculo de los coeficientes:

vcoef = zeros(Nx,Ny,Nx,Ny);
for i = 0:Nx-1
    progress_bar(i + 1, Nx, 40, 'Computing coefficients');
    for j = 0:Ny-1
        for k = 0:Nx-1
            for l = 0:Ny-1
                vcoef(i+1,j+1,k+1,l+1) = (1/(Lx*Ly)) * ...
                    integral_ijkl(v1grid, v2grid, i, j, k, l, Nx, Ny, Lx, Ly);
            end
        end
    end
end

%% Almacenamiento de los datos:

N2 = Nx*Ny;

A = zeros(N2, N2);
A_row = zeros(N2^2, 1);

B = zeros(N2,N2);
B_row = zeros(N2^2, 1);

Lambda = zeros(N2, N2);
Lambda_row = zeros(N2^2, 1);

u0_proj_row = zeros(N2,1);

r = 0;
for i = 1:Nx
    progress_bar(i, Nx, 40, 'Storing Vectorized Forms');
  for j = 1:Ny
    m = i+(j-1)*Nx;
    u0_proj_row(m,1) = u0_proj(i,j);
    for k = 1:Nx
      for l = 1:Ny
        n = k+(l-1)*Nx; %Corrected: k+(l-1)*Ny
        r = r+1; %Va contando de acuerdo a 
        %Matrix A:
        A(m,n) = vcoef(i,j,k,l);
        A_row(r,1) = A(m,n);
        %Matrix B & Lambda:
        if i == k && j == l
            lambda_ij = (pi^2) * (i^2) / (Lx^2) + (pi^2) * (j^2) / (Ly^2);
            B(m, n) = lambda_ij^(-gamma);
            B_row(r, 1) = B(m,n);
            %
            Lambda(m,n) = lambda_ij;
            Lambda_row(r, 1) = Lambda(m,n);
        end
      end
    end
  end
end

%Matrix A
figure(1)
contour(A,100); colorbar;
set(gca,'Ydir','reverse')

print(gcf, 'Data/MatrizAGrafica','-depsc')
save('Data/MatrizA.mat','A','Lx','Ly','Nx','Ny')
csvwrite('Data/MatrixA.csv',A_row)
save Data/MatrixA.dat A_row -ascii

%Matrix B
figure(2)
contour(B,100); colorbar
set(gca,'Ydir','reverse')
print(2, 'Data/MatrizBGrafica','-depsc')
save('Data/MatrizB.mat','B','Lx','Ly','Nx','Ny')
csvwrite('Data/MatrixB.csv',B_row)
save Data/MatrixB.dat B_row -ascii

%Matrix Lambda
figure(3)
contour(Lambda,100); colorbar;
set(gca,'Ydir','reverse')
print(3, 'Data/MatrizLambdaGrafica','-depsc')
save('Data/MatrizLambda.mat','Lambda','Lx','Ly','Nx','Ny')
csvwrite('Data/MatrixLambda.csv',Lambda_row)
save Data/MatrixLambda.dat Lambda_row -ascii

%Initial condition projected
save('Data/u0_proj_row.mat','u0_proj_row','Nx','Ny')
csvwrite('Data/u0_proj_row.csv',u0_proj_row)
save Data/u0_proj_row.dat u0_proj_row -ascii

u0_proj_row;
