function Evolution

%% Parameters:

Nx = 50; Ny = 50; %Numeros de puntos en la malla
Lx = 5; Ly = 5;

beta = 0.1;
theta = 1;
T = .5;
dt = 10^(-5); %Tamanio de paso
%T = 10000*dt;

%% Load projected initial condition

load('Data/u0_proj_row.dat')
%Loads u0_proj_row
%Nx and Ny

load('Data/MatrixLambda.dat')
%Loads Lambda_row

% load('MatrixB.dat')
% %Loads B_row

load('Data/MatrixA.dat')
%Loads A_row

N2 = Nx*Ny;
%vcoef = zeros(Nx,Ny,Nx,Ny);
%u0_proj = zeros(Nx,Ny);

A = zeros(N2,N2);
%B = zeros(N2,N2);
Lambda = zeros(N2,N2);

r = 0;
for i = 1:Nx
    for j = 1:Ny
        m = i+(j-1)*Nx; %Corrected: m = i+(j-1)*Ny
        %u0_proj(i,j) = u0_proj_row(m,1);
        for k = 1:Nx
            for l = 1:Ny
                n = k+(l-1)*Nx; %Corrected: k+(l-1)*Ny
                r = r+1; %Va contando de acuerdo a 
                %Matrix A:
                %vcoef(i,j,k,l) = A_row(r,1);
                A(m,n) = MatrixA(r,1);
                %Matrix B & Lambda:
                if i == k && j == l
                    %B(m,n) = MatrixB(r,1);
                    %
                    Lambda(m,n) = MatrixLambda(r,1);
                end
            end
        end
    end
end

%% Time evolution

u0_proj_row_np1 = u0_proj_row;
u0_proj_row_n = u0_proj_row;

M = beta*Lambda+theta*A;
t = 0;
while t < T
    t
    u0_proj_row_np1 = u0_proj_row_n-M*u0_proj_row_n*dt;

    u0_proj_row_n = u0_proj_row_np1;

    t = t+dt;
end

%Pass row to matrix form:

u0_proj_np1 = zeros(Nx,Ny);
for i = 1:Nx
    for j = 1:Ny
        m = i+(j-1)*Nx; %Corrected: m = i+(j-1)*Ny
        u0_proj_np1(i,j) = u0_proj_row_np1(m,1);
    end
end

u0_comp_np1_grid = comp(u0_proj_np1,Nx,Ny,Lx,Ly);

plot(log(abs(u0_proj_np1(:,15))))
stop

%% Plots:


xgrid = zeros(Nx,Ny);
ygrid = zeros(Nx,Ny);
v1grid = zeros(Nx,Ny);
v2grid = zeros(Nx,Ny);
u0grid = zeros(Nx,Nx);
for iy = 1:Ny
    for ix = 1:Nx
        xgrid(ix,iy) = Lx/Nx*(ix-1/2);
        ygrid(ix,iy) = Ly/Ny*(iy-1/2);
        v1grid(ix,iy) = v1(xgrid(ix,iy),ygrid(ix,iy));
        v2grid(ix,iy) = v2(xgrid(ix,iy),ygrid(ix,iy));
        
        xi = xgrid(ix,iy);
        yj = ygrid(ix,iy);
        u0grid(ix,iy) = u0(xi,yj,Lx,Ly);
    end
end

figure(1)
subplot(2,2,1)
hold on
title('u_0')
quiver(xgrid,ygrid,v1grid,v2grid)
contour(xgrid,ygrid,u0grid,100); colorbar;
subplot(2,2,2)
hold on
title('u at t = .5')
quiver(xgrid,ygrid,v1grid,v2grid)
contour(xgrid,ygrid,u0_comp_np1_grid,100); colorbar;
axis([0,Lx,0,Ly])

