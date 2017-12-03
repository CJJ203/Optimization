clear;
close all;
warning off;

% Example 1:
syms x1 x2 x3;
f = (x1-4)^4 + (x2-3)^2 +4*(x3+5)^4;
[x,fx] = SteepestDescent(f, [x1,x2,x3], [4,2,-1], 0.0001);

% Example 2:
syms x1 x2;
f = 3/2*x1^2+1/2*x2^2-x1*x2-2*x1;
var = [x1 x2];
x0 = [0 0];
e = 1e-2;
[x fx] = SteepestDescent(f, var, x0, e);