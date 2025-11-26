% Clear workspace and command window for clean execution
clear;
close all;
clc;

% Parameters definition
Lx = 5;
Ly = 5;
gamma = 2;
Nx = 50;
Ny = 50;

% Initialize grid arrays using vectorization where possible
[xgrid, ygrid] = ...
    meshgrid( ...
        linspace( ...
            Lx / (2 * Nx), ...
            Lx - Lx/(2 * Nx), Nx ...
        ), ...
        linspace( ...
            Ly / (2 * Ny), ...
            Ly - Ly/(2 * Ny), Ny ...
        ) ... 
    );
xgrid = xgrid';
ygrid = ygrid';

% Evaluate velocity and initial condition on the grid (vectorized)
v1grid = v1(xgrid, ygrid, Lx, Ly);
v2grid = v2(xgrid, ygrid, Lx, Ly);
u0grid = u0(xgrid, ygrid, Lx, Ly);


%% Initial condition projection
u0_proj = compute_initial_projection(u0grid, Nx, Ny, Lx, Ly);

u0_comp_grid = comp(u0_proj,Nx,Ny,Lx,Ly);

figure(1)
subplot(2, 2, 1)
hold on
quiver(xgrid, ygrid, v1grid, v2grid)
contour(xgrid, ygrid, u0grid, 100); 
colorbar;
subplot(2, 2, 2)
hold on
quiver(xgrid, ygrid, v1grid, v2grid)
contour(xgrid, ygrid, u0_comp_grid, 100); 
colorbar;
axis([0,Lx,0,Ly])
print('CampoDeVectores','-depsc',figure(1))
clf 

%stop

%% Compute velocity coefficients
vcoef = compute_velocity_coefficients(v1grid, v2grid, Nx, Ny, Lx, Ly);

%% Build and store matrices
[A, B, Lambda, u0_proj_row] = ...
    build_matrices(vcoef, u0_proj, Nx, Ny, Lx, Ly, gamma);

%% Visualizations and data storage
visualize_and_save_results(A, B, Lambda, u0_proj_row, Lx, Ly, Nx, Ny);


