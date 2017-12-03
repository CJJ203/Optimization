%%%%%%%%%%%%%%% Steepest Descent Method %%%%%%%%%%%%%%%
% Author:   ²Ì¼Ñ¼Ñ
% ID:       SZ1716003
% Date:     2017-12-01
%
% Input:
%   f:      objective function
%   var:    row vector, variables in function f
%   x0:     row vector, start point
%   eps:    accuracy
%
% Output:
%   x:      the local minimal point
%   minf:   the local minimal value of f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, minf] = SteepestDescent(f, var, x0, eps)

if nargin == 3
    eps = 1.0e-3;
end

tolerance = 1;

iter = 1;

message = @(iter,x,fx,ss) fprintf('iter = %3d: x = %-24s, f(x) = %12f, step_size = %f\n', iter, mat2str(x,6), fx, ss);


while tolerance > eps
    
    f_x = double(subs(f, var, x0));
    
    gradf = jacobian(f, var);
    grad_x0 = double(subs(gradf, var, x0));
    tolerance = norm(grad_x0);
    
    syms a;
    x_temp = x0 - a * grad_x0;
    f_temp = subs(f, var, x_temp);
    h = diff(f_temp, a);
    
    if h ~= 0,
        opt_alpha = solve(h, 'Real', true);
        opt_alpha = double(opt_alpha(1));
    else
        break;
    end
    
    message(iter, x0', f_x, double(opt_alpha));
    
    x1 = x0 - opt_alpha * grad_x0;
    x0 = double(x1);
    iter = iter + 1;
end

x = x0;
minf = double(subs(f,var,x));
fprintf('Answer: x = %-24s, f(x) = %f\n\n', mat2str(x0'), minf);
end
