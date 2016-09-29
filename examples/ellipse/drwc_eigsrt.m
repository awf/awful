function [Vout,Dout] = drwc_eigsrt(M)
% EIGSRT     Sorted eigensystem
%            EIGSRT is a wrapper for EIG which sorts the eigensystem
%            by ascending eigenvalue.
 
[V, D] = eig(M);
[D, Di] = sort(diag(-D));
if D(2) < 0
  D=-D;
else
  D=flipud(D);
  Di=flipud(Di);
end
 
 
if nargout <= 1
  Vout = D;
else
  Vout = V(:,Di);
  Dout = D;
end
