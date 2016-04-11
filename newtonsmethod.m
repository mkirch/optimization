function [dk,xkplusone] = newtonsmethod(f,x0,max)
% f - function
% x0 - initial vector starting point
xi=x0;

a=1;
k=0;

grad = gradient(f);
hess = hessian(f);

while k<max
    fxi = subs(f,transpose(symvar(f)),xi);
    double(fxi)
    gxi = subs(grad,transpose(symvar(f)),xi);
    hxi = subs(hess,transpose(symvar(f)),xi);
    dk = -inv(hxi)*gxi;
    xplusone = xi+dk*a;
    xi=xplusone;
    k=k+1;
end
end

