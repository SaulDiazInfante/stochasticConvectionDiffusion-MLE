function CoefficientsUijkl

% load('Coeficientes.mat')
% vcoef
% stop

%% Definicion de campo de velocidades:

L1 = 5; L2 = 5;

% hx = @(x,y) sin(x).*sin(x);
% integral2(hx,0,2*pi,0,2*pi)

%% Calculo de los coeficientes:

Nx = 10;
Ny = 10;

vcoef = zeros(Nx,Ny,Nx,Ny);

v1xy =@(x,y) v1(x,y);
v2xy =@(x,y) v2(x,y);
for i = 1:Nx
    i
    hi =@(x,y) cos(pi*i/L1*x);
    hi_prime =@(x,y) -pi*i/L1*sin(pi*i/L1*x);
    for j = 1:Ny
        hj =@(x,y) cos(pi*j/L1*y);
        hj_prime =@(x,y) -pi*j/L1*sin(pi*j/L1*y);
        for k = 1:Nx
            hk =@(x,y) cos(pi*k/L1*x);
            for l = 1:Ny
                hl =@(x,y) cos(pi*l/L1*y);

                integrand_xy =@(x,y) v1xy(x,y).*hi_prime(x,y).*hj(x,y).*hk(x,y).*hl(x,y)+v2xy(x,y).*hi(x,y).*hj_prime(x,y).*hk(x,y).*hl(x,y);

                %integrand_xy =@(x,y) v1xy(x,y).*(-pi*i/L1*sin(pi*i/L1*x)).*(cos(pi*j/L1*y)).*(cos(pi*k/L1*x)).*(cos(pi*l/L1*y))+v2xy(x,y).*(cos(pi*i/L1*x)).*(-pi*j/L1*sin(pi*j/L1*y)).*(cos(pi*k/L1*x)).*(cos(pi*l/L1*y));

                vcoef(i,j,k,l) = integral2(integrand_xy,0,L1,0,L2);

            end
        end
    end
end

N2 = Nx*Ny;

A = zeros(N2,N2);
Arow = zeros(N2^2,1);

r = 0;
for i = 1:Nx
    for j = 1:Ny
        m = i+(j-1)*Ny;
        for k = 1:Nx
            for l = 1:Ny
                n = k+(l-1)*Ny;
                A(m,n) = vcoef(i,j,k,l);
                r = r+1; %Va contando de acuerdo a 
                 Arow(r,1) = A(m,n);
            end
        end
    end
end

figure(1)
%contour(vcoef(:,:,Nk,Nl),100); colorbar
contour(A,100); colorbar
%surf(vcoef)
print('MatrizAGrafica','-depsc',figure(1))
save('MatrizA.mat','A','L1','L2','Nx','Ny')
csvwrite('MatrixA.csv',Arow)
save MatrixA.dat Arow -ascii
stop

%% Campo de velocidades

N1 = 30;
N2 = 30;

xgrid = zeros(N1,N2);
ygrid = zeros(N1,N2);
v1grid = zeros(N1,N2);
v2grid = zeros(N1,N2);

for iy = 1:N2
    for ix = 1:N1
        xgrid(ix,iy) = L1/N1*(ix-1/2);
        ygrid(ix,iy) = L2/N2*(iy-1/2);
        v1grid(ix,iy) = v1(xgrid(ix,iy),ygrid(ix,iy));
        v2grid(ix,iy) = v2(xgrid(ix,iy),ygrid(ix,iy));
    end
end

figure(1)
quiver(xgrid,ygrid,v1grid,v2grid)
axis([0,L1,0,L2])

print('CampoDeVectores','-depsc',figure(1))

