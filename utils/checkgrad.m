function err = checkgrad(f,df)
% CHECKGRAD   verify an function in R3 has correct analytic gradient
%
% err = checkgrad(f,df) returns error between finite-differencing approx to
%  grad f and the function df. f is a function handle taking a stack of m 3-cpt
%  col vecs and returning a row vec of m values. df only need map a single
%  col vec to a col vec gradient vector.
%
% a self-test is done with no inputs

% Barnett 7/18/16

if nargin==0, test_checkgrad; return; end
x0 = [0.4;0.6;-0.3];
eps = 1e-5;
x = repmat(x0,[1,7]);  % set up coords to eval f at
for i=1:3, v = zeros(3,1); v(i) = 1;  % unit vec
  x(:,2*i) = x(:,2*i)-eps*v; x(:,2*i+1) = x(:,2*i+1)+eps*v;
end
u = f(x);
du = [u(3)-u(2);u(5)-u(4);u(7)-u(6)]/(2*eps);  % hmm, I guess u(1) ignored
err = norm(du-df(x0));

%%%%%%%%
function test_checkgrad
a = [1;-1.6;0.7];          % the const grad
f = @(x) 7.0 + a(1)*x(1,:) + a(2)*x(2,:) + a(3)*x(3,:);  % linear func
df = @(x) repmat(a,[1 size(x,2)]);
checkgrad(f,df)
