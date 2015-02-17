function J = au_jacobian_fd(f,x,varargin)
% AU_JACOBIAN_FD Finite-difference jacobian
%     J = AU_JACOBIAN_FD(f,x, opts)
%     With opts an au_opts structure defaulting to
%         delta = 1e-4      -- Added to each element of x
%         verbose = 1       -- Be verbose
%         

% awf, jun13

if nargin == 0
    %% Test case
    f = @(x) [norm(x); norm(x).^3];
    jac = @(x1,x2) ...
        [   x1/(x1^2 + x2^2)^(1/2),   x2/(x1^2 + x2^2)^(1/2)
        3*x1*(x1^2 + x2^2)^(1/2), 3*x2*(x1^2 + x2^2)^(1/2)];
    x = rand(2,1);
    fdJ = au_jacobian_fd(f, x);
    gtJ = jac(x(1),x(2));
    au_test_begin
    au_test_equal gtJ fdJ 1e-4
    au_test_end
    return
end

opts = au_opts('delta=1e-4;verbose=1', varargin{:});

fx = f(x);
n = length(fx);
p = length(x);
J = zeros(n,p);
au_assert_equal p numel(x)
for k=1:p
    e = zeros(size(x));
    e(k) = opts.delta;
    scale = 1/((x(k)+e(k))-(x(k)-e(k)));
    J(:,k) = (f(x+e) - f(x-e))*scale;
end
