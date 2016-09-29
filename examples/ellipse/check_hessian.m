function check_hessian(f, x)

%% test
if nargin == 0
  test
  return
end

%% fun
delta = 1e-7;

fprintf('check hessian and gradient:');
Jfd = fdgrad(f,x,delta);
Hfd = fdgrad(@(x) grad(f, x), x, delta);

[~,J,H] = f(x);

H = full(H);
Jerr = norm(J - Jfd, 'fro');
fprintf('J err = %.06f/%.1f, ', Jerr, norm(J));

Herr = norm(H - Hfd, 'fro');
fprintf('H err = %.06f/%.1f = %.06f\n', Herr, norm(H), Herr/norm(H));

if 1
if Jerr > 1e-4
  disp('J|Jfd|Jdiff = ');
  au_prmat(J,Jfd,J-Jfd)
end

if Herr > 1e-5 * norm(H(:))
  disp('H|Hfd|Hdiff = ');
  au_prmat(H, Hfd, H-Hfd);
end
end

function g = grad(f, x)
[~,g] = f(x);

%%
function g = fdgrad(f, x, delta)
x = x(:);
fx = f(x);
d = size(x,1);
n = size(fx,1);
g = zeros(d,n);
for i=1:d
  xp = x; xp(i) = xp(i) + delta;
  xn = x; xn(i) = xn(i) - delta;
  g(i,:) = (f(xp) - f(xn))/(xp(i) - xn(i));
end

%% 
function test
disp('** test check_hessian **');
check_hessian(@testfun, [.1 .3]);
check_hessian(@testfun, [1 1]);
check_hessian(@testfun, randn(1,2));


%%
function [e, J, H] = testfun(p)
x = p(1);
y = p(2);
e =(10*y-10*x^2)^2+(1-x)^2;
J =[
 -40*(10*y-10*x^2)*x-2+2*x
             200*y-200*x^2
];
H = [
 1200*x^2-400*y+2,           -400*x
           -400*x,              200
];
