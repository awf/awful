function nll = example_gmm(alphas, mus, Ls, x)
% EXAMPLE_GMM  Evaluate GMM negative log likelihood for a single point x
%      Inverse covariance lower triangular roots are in dxdxK array Ls.
%      That means kth covariance matrix is inv(Ls(:,:,k)*Ls(:,:,k)')
%      This will be mexed by au_ccode, so doesn't need to be super fast.
%      Answer is scaled by (2*pi)^d/2


% weight[k] = exp(alphas[k])/sum(exp(alphas))
% mahal[k] = Ls[k]*(mus[k] - x)
% log(sum_k weight[k] * det(Ls[k]) * exp(-0.5*sumsq(mahal[k])))
% =log(sum_k exp(alphas[k])/sum(exp(alphas)) * exp(log(det(Ls[k])) * exp(-0.5*sumsq(mahal_k)))
% =log(1/sum(exp(alphas)) * sum_k { exp(alphas[k]) * exp(log(det(Ls[k])) * exp(-0.5*sumsq(mahal_k)))
% =log(1/sum(exp(alphas)) * sum_k exp(alphas[k]) + log(det(Ls[k])) - 0.5*sumsq(mahal_k)))
% =-log(sum(exp(alphas)) + 
%       log(sum_k exp(alphas[k] + log(det(Ls[k])) - 0.5*sumsq(mahal_k)))
% =-log(sum(exp(alphas)) + 
%       log(sum_k exp(alphas[k] + log(prod(diag(Ls[k]))) - 0.5*sumsq(mahal_k)))
% =-log(sum(exp(alphas)) + 
%       log(sum_k exp(alphas[k] + log(prod(exp(diag(Ls0[k])))) - 0.5*sumsq(mahal_k)))
% =-log(sum(exp(alphas)) + 
%       log(sum_k exp(alphas[k] + sum(diag(Ls0[k])) - 0.5*sumsq(mahal_k)))

if nargin == 0
    %% test
    K  = 3;
    d = 7;
    alphas = randn(K,1);
    mus = randn(d,K);
    Ls = randn(d,d,K);
    x = randn(d,1);
    
    weights = exp(alphas)/sum(exp(alphas));
    likes = zeros(K,1);
    likes1 = zeros(K,1);
    for k=1:K
        L0 = Ls(:,:,k);
        L = tril(L0, -1) + diag(exp(diag(L0)));
        iCov = L'*L;
        dx = x - mus(:,k);
        likes1(k) = weights(k)*det(iCov)^0.5*exp(-0.5*dx'*iCov*dx);

        logw = alphas(k)-au_logsumexp(alphas);
        mahal = L*dx;
        likes(k) = logw + log(prod(exp(diag(L0)))) - 0.5*(mahal'*mahal);
        likes(k) = logw + sum(diag(L0)) - 0.5*(mahal'*mahal);
    end
    au_test_equal log(likes1) likes 1e-7
    nll0 = au_logsumexp(log(likes1));
    
    nll1 = example_gmm(alphas, mus, Ls, x);

    au_test_equal nll0 nll1 1e-7
    return
end

K = length(alphas);
lse = zeros(K,1, 'like', alphas);
for k=1:K
  L0 = Ls(:,:,k);
  L = tril(L0, -1) + diag(exp(diag(L0)));
  
  mahal = L*(mus(:,k) - x);
  lse(k) = alphas(k) + sum(diag(L0)) - 0.5*(mahal'*mahal);
end
if isa(lse, 'sym')
    nll = au_logsumexp(lse) - au_logsumexp(alphas);
else
    logsumexp = @(x) log(sum(exp(x)));
    nll = logsumexp(lse) - logsumexp(alphas);
end
