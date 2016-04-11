function [x,fval,g,nfe,nge,xs] = ...
    steepest_ls(func,x,alpha,c1,rho,tol,itmax,trace)
% [x,fval,gval,nfe,nge,xs] = ...
%	steepest_ls(func,x,alpha,c1,rho,tol,itmax,trace)
%
% Performs steepest descent search from x with
% initial step-length parameter alpha.
% Line search termination criteria are the sufficient decrease
% condition:
%   f(x+alpha*p) <= f(x) + c1*alpha*g'*p
% where p = -g, the steepest descent direction
% Backtracking with parameter RHO (between 0 and 1) is used as
% the line-search method.
%    func   -- function name; 
%              for [fval,gval] = func(x), fval is the function value at x,
%              and gval is the gradient value
%    tol    -- stopping tolerance for gradient
%    itmax  -- maximum number of iterations
% Returns x and func(x), # function eval'ns (nfe) and # gradient
% eval'ns (nge), and also the gradient at x (gval).
% If trace is non-zero, then information will be printed about the process

% Initialization
if nargin <= 7
  trace = 0
end
nfe = 0;
nge = 0;
alpha0 = alpha;

xs = x;

for k = 1:itmax
  % Choose steepest descent direction.
  [fval0,g] = feval(func,x);
  if trace > 0
    k
    x
    g
  end
  nge = nge + 1;
  if norm(g) < tol
      fval = fval0;
      break
  end
  [m,n] = size(g);
  if m < n
    g = g';
  end
  p = -g;
  slope = g'*p;

  alpha = alpha0;
  % Check sufficient decrease criterion
  % fval0 = feval(func,x);
  fval1 = feval(func,x+alpha*p);
  nfe = nfe + 1;
  
  % Backtrack while criterion not satisfied
  while fval1 > fval0 + c1*alpha*slope
    alpha = rho*alpha;
    fval1 = feval(func,x+alpha*p);
    if trace > 1
      alpha
      fval1
    end
    nfe = nfe + 1;
  end

  x = x+alpha*p;
  if ( trace > 0 )
    xs = [xs, x];
  end
  fval = fval1;
  fval0 = fval;
  if trace > 0
    fprintf('Steepest descent step with alpha = %g\n',alpha);
  end
end
