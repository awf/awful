function [x, lambda] = mineig(A)
% MINEIG     Return eigenvector of A corresponding to smallest eigenvalue 

[v, evalues] = eig(A);
evalues = diag(evalues);
[lambda,dminindex] = min(evalues);
x = v(:,dminindex(1));
