function x = cmin(A,B)
% CMIN       Constrained minimization.
%       Minimize x' A x subject to x'B x < 0
%       A non-negative definite, B real symmetric
[n, m] = size(A);
[U, Da] = eig(A); da = diag(Da);

[V, Db] = eig(B);
db = diag(Db);
[Ddum, Di] = sort(-abs(db));
db = db(Di);
Db = diag(db);
V = V(:,Di);

debug = 1;
if debug
  if any(da < -sqrt(eps))
    A
    EigenVectors_A = U
    EigenValues_A = da
    error('A not non-negative definite');
  end
  if any(imag(db))
    error('B not real symmetric');
  end
end 

rank_a = sum(da > eps);
%if rank_a < n
%  A
%  EigenVectors_A = U
%  EigenValues_A = da
%  error('rank(A) < n');
%end

rank_b = sum(abs(db) > eps);
if rank_b == n
  invB = V * diag(1./db) * V';
  x = mineig(invB * A);
else

  M = V'*U*Da*U'*V;
  I1 = 1:rank_b; I2 = rank_b+1:n;
  M11 = M(I1, I1); M12 = M(I1, I2);
  M21 = M(I2, I1); M22 = M(I2, I2);
  invDb11 = diag(1./db(I1));
  
  invM22_M21 = inv(M22) * M21;
  N = invDb11 * (M11 - M12 * invM22_M21);
  y1 = mineig(N);
  y2 = -invM22_M21 * y1;
  y = [y1;y2];
  x = V*y;

  Constraint = x'*B*x;
  if Constraint > 0
    % keyboard
    x = inf*x;
  end

end
