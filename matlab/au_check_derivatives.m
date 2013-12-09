function emax = au_check_derivatives(f,x,J,delta, tol)
% AU_CHECK_DERIVATIVES Finite-difference derivative check
%     emax = au_check_derivatives(f,x,J[,delta[,tol]])
%     Default delta = 1e-4, tol = 1e-7

if nargin == 0
    %% Test case
    f = @(x) [norm(x); norm(x).^3];
    jac = @(x1,x2) ...
        [   x1/(x1^2 + x2^2)^(1/2),   x2/(x1^2 + x2^2)^(1/2)
        3*x1*(x1^2 + x2^2)^(1/2), 3*x2*(x1^2 + x2^2)^(1/2)];
    disp('Should be silent, returning emax:');
    au_check_derivatives(f, [.2 .3]', jac(.2, .3))
    disp('Should throw, jac evaluated at wrong place: ');
    au_check_derivatives(f, [.2 .3]', jac(.2, .31))
    return
end

if nargin < 4
    delta=1e-4;
end

if nargin < 5
    tol=1e-7;
end

scale = 1/2/delta;
[~,p] = size(J);
au_assert_equal p numel(x)
% Check derivatives OK
fdJ=0*J;
for k=1:p
  e = zeros(size(x));
  e(k) = delta; 
  fdJ(:,k) = (f(x+e) - f(x-e))*scale;
end
emax = max((fdJ(:) - J(:)).^2);
au_assert_equal fdJ J 1e-7
