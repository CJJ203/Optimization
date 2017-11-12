function [x, minf] = SteepestDescent(f, x0, var, eps)

format long;
if nargin == 3
    eps = 1.0e-6;
end

syms l;
tol = 1;

while tol > eps
    gradf = -jacobian(f, var);
    v = fnval(gradf, var, x0);
    tol = norm(v);
    y = x0 + l * v;
    yf = fnval(f, var, y);
%     [a,b] = minJT(yf,0,0.1);
    xm = GoldenSearch(yf,0,0.1);
    x1 = x0 + xm*v;
    x0 = x1;
end

x = x1;
minf = fnval(f,var,x);
format short;

