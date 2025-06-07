function u0xy = u0(x,y,Lx,Ly)

frac = 1/4 ;

rad = sqrt( ((x-Lx/4)/(frac*Lx))^2+((y-0.6*Ly)/(frac*Ly))^2 );

u0xy = 0;
if rad <= 1
    u0xy = (1+cos(pi*rad))/2;
end



