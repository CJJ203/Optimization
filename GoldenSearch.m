function [x,minf] = GoldenSearch(f,a,b,eps)

if nargin == 3
    eps = 1.0e-6; % default accuracy
end

l = a + 0.382*(b-a);
u = a + 0.618*(b-a);
k = 1;
tol = b-a;

while tol>eps && k<1000
    fl = subs(f, findsym(f), l);
    fu = subs(f, findsym(f), u);
    if fl > fu
        a = l;
        l = u;
        u = a + 0.618*(b-a);
    else
        b = u;
        u = l;
        l = a + 0.382*(b-a);
    end
    k = k+1;
    tol = abs(b-a);
end

if k == 1000
    disp('Cannot find the minimizer!');
    x = NaN;
    minf = NaN;
    return;
end

x = (a+b)/2;
minf = double(subs(f,x));

