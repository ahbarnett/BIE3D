% tester for MEX extrap routine. Barnett 12/23/16

%f = @(x) x.^2;
%f = @(x) sin(0.3*x);
f = @(x) exp(0.1*x);
m = 10;
w = extrap(m);
w*f(0:m)' - f(m+1)
