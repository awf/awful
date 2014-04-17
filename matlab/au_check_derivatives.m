function emax = au_check_derivatives(f,x,J,delta, tol, timeout, verbose)
% AU_CHECK_DERIVATIVES Finite-difference derivative check
%     emax = au_check_derivatives(f,x,J[,delta[,tol]])
%     Default delta = 1e-4, tol = 1e-7

if nargin == 0
    %% Test case
    f = @(x) [norm(x); norm(x).^3];
    jac = @(x1,x2) ...
        [   x1/(x1^2 + x2^2)^(1/2),   x2/(x1^2 + x2^2)^(1/2)
        3*x1*(x1^2 + x2^2)^(1/2), 3*x2*(x1^2 + x2^2)^(1/2)];
    disp('Should be silent...');
    au_check_derivatives(f, [.2 .3]', jac(.2, .3));
    disp('was it?');
    au_test_should_fail au_check_derivatives(f,[.2,.3]',jac(.2, .31))
    return
end

if nargin < 4
    delta=1e-4;
end

if nargin < 5
    tol=1e-7;
end

if nargin < 6
    timeout = inf;
end

if nargin < 7
    verbose = 1;
end


scale = 1/2/delta;
[~,p] = size(J);
if verbose
    fprintf('au_check_derivatives: dim %d, time at most %.1f seconds...', p, timeout);
end
au_assert_equal p numel(x)
% Check derivatives OK
fdJ=0*J;
t = clock;
emax = 0;
ks = randperm(p);
for ki=1:p
    k = ks(ki);
    e = zeros(size(x));
    e(k) = delta;
    fdJ(:,k) = (f(x+e) - f(x-e))*scale;
    err = max(abs(fdJ(:,k) - J(:,k)));
    if err > tol
        err = full(err);
        error('awful:check_derivatives', 'au_check_derivatives: Error on parameter %d = %g', k, err);
    end
    emax = max(emax, err);
    if etime(clock, t) > timeout
        if verbose
            fprintf('timeout, checked %d]', ki);
        end
        break
    end
end
%au_assert_equal('fdJ', 'J', tol)
