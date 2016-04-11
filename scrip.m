syms x y;

%Function and Inital X_0
f = symfun(y^2-(x^2)*y - 3*y+x^4 + -x^3+2*x*y-x,[x y]);
x0=[0;1];

newtonsmethod(f,x0,5)

