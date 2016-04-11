function [xi,fxi,gxi] = newtonsmethod(f,x0,max,c1,tol)
% f - function
% x0 - initial vector starting point
% max - maximum number of iterations
% c1 - Armijo constant
format long
xi=x0;
k=0;
grad = gradient(f);
hess = hessian(f);

while k<max
    fxi = subs(f,transpose(symvar(f)),xi);
    gxi = subs(grad,transpose(symvar(f)),xi);
    hxi = subs(hess,transpose(symvar(f)),xi);
    % If hess positive definite, use Newton's. 
    % If hess not positive definite, use Steepest Descent.
    [R,p] = chol(hxi);
    if p==0
        dk = -inv(hxi)*gxi;
    else 
        dk = -gxi;
    end
    
    % Determine Step Length Using Armijo-Backtracking Algorithm
    a = 1;
    while (subs(f,transpose(symvar(f)),xi+a*dk)) > (fxi+c1*a*transpose(dk)*gxi)
        a = a/2;
    end
    
    % Next x
    xk1 = xi+dk*a;
    
    % Stop or update
    if abs(xi - xk1) < tol
        xi = double(xk1);
        fxi = double(subs(f,transpose(symvar(f)),xi));
        gxi = double(subs(grad,transpose(symvar(f)),xi));
        return            
    else
        xi = xk1; 
        k=k+1;
    end
end
end

