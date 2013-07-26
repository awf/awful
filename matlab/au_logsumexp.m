function [l, Jacobian] = au_logsumexp(M)
% au_LOGSUMEXP   Compute log(sum(exp(M))) stably
%                   au_logsumexp(M) = log(sum(exp(M))) 
%                 but avoids under/overflow.

% awf, may13

if nargin == 0
  %% test
  Ms = {
    -rand(7,3)
    -rand(1,7)
    -rand(7,1)
    randn(7,3,4)
    [0 0 0 -5 -10 -15 -20 -10 -5]-1000
    [0 0 0 -5 -10 -15 -20 -10 -5]+900
    };
  for k = 1:length(Ms)
    M = Ms{k};
    fprintf('test %d size', k); fprintf(' %d', size(M)); fprintf('\n');
    lseM = log(sum(exp(M)));
    au_lseM = au_logsumexp(M);
    q = @(x) squeeze(x);
    au_prmat(q(lseM), q(au_lseM));
  end
  return
end

%% Main body
sz = size(M);
if prod(sz) == max(sz)
  A = max(M);
  l = log(sum(exp(M - A))) + A;
else
  A = max(M);
  l = bsxfun(@plus, log(sum(exp(bsxfun(@minus, M, A)))), A);
end
