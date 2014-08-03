function emax = au_check_derivatives(f,x,J,varargin)
% AU_CHECK_DERIVATIVES Finite-difference derivative check
%     emax = au_check_derivatives(f,x,J, opts)
%     With opts an au_opts structure defaulting to
%         delta = 1e-4      -- Added to each element of x
%         tol = 1e-7        -- Want derivatives accurate to this tolerance
%         timeout = inf     -- Spend at most "timeout" seconds checking
%         verbose = 1       -- Be verbose
%         PatternOnly = 0   -- Used for checking JacobPattern
%         

% awf, jun13

if nargin == 0
    %% Test case
    f = @(x) [norm(x); norm(x).^3];
    jac = @(x1,x2) ...
        [   x1/(x1^2 + x2^2)^(1/2),   x2/(x1^2 + x2^2)^(1/2)
        3*x1*(x1^2 + x2^2)^(1/2), 3*x2*(x1^2 + x2^2)^(1/2)];
    au_check_derivatives(f, [.2 .3]', jac(.2, .3));
    fprintf('Should be silent [...');
    au_check_derivatives(f, [.2 .3]', jac(.2, .3), 1e-5, 1e-5, inf, 0);
    disp('] was it?');
    au_test_should_fail au_check_derivatives(f,[.2,.3]',jac(.2, .31))
    return
end

opts = au_opts('delta=1e-4;tol=1e-7;timeout=60;verbose=1;PatternOnly=0', varargin{:});

[~,p] = size(J);
if opts.verbose
    fprintf('au_check_derivatives: dim %d, time at most %.1f seconds...', p, opts.timeout);
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
    e(k) = opts.delta;
    scale = 1/((x(k)+e(k))-(x(k)-e(k)));
    fdJ(:,k) = (f(x+e) - f(x-e))*scale;
    % Used for checking JacobPattern
    if opts.PatternOnly
        fdJ(:,k) = fdJ(:,k) ~= 0;
        err = max(fdJ(:,k) - J(:,k));
    else
        err = max(abs(fdJ(:,k) - J(:,k)));
    end
    if err > opts.tol
        err = full(err);
        error('awful:check_derivatives', 'au_check_derivatives: Error on parameter %d = %g', k, err);
    end
    emax = max(emax, err);
    if etime(clock, t) > opts.timeout
        if opts.verbose
            fprintf('timeout, checked %d/%d]\n', ki, p);
        end
        break
    end
end
if opts.verbose
    if opts.PatternOnly
        fprintf(' %d/%d extra zeros in Pattern ', full(sum(J(:)-fdJ(:))), nnz(J));
    end
    fprintf('all OK\n');
end
%au_assert_equal('fdJ', 'J', tol)
