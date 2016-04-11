syms x y;

%Function and Initial X_0
f = symfun(y^2-(x^2)*y - 3*y+x^4 + -x^3+2*x*y-x,[x y]);
x0=[2;2];

[xi,fxi,gxi] = newtonsmethod(f,x0,50,10^-4,0.00000000001)

