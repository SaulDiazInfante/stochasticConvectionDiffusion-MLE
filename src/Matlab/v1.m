function v1xy = v1(x,y)

r =1 ; L1 = 5; L2 = 5;
v1xy = ( y+x+2*cos(x/2)+sin(x/2) ).*exp(-.6*((x-L1/2).^2+(y-L2/2).^2).^r);