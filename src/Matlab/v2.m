function v2xy = v2(x,y)

r =1 ; L1 = 5; L2 = 5;
v2xy = ( (-8/17)*exp(x/2).*cos(2*x)+(2/17)*exp(x/2).*sin(2*x)-y ).*exp(-.6*((x-L1/2).^2+(y-L2/2).^2).^r);