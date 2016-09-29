function [r, J] = au_grad_roots(coefs)
% AU_GRAD_ROOTS  Find polynomial roots and gradients wrt params

if nargin == 0
  %% Test
  coef_dups = [1 -8 23 -28 12];
  coef_close = [1.0000   -8.1000   23.6000  -29.1000   12.6000];
  coef = coef_close
  clf
  t = 0.5:.01:3.5;
  plot(t, polyval(coef, t));
  hold on
  plot(t,0*t, 'r-');
  axis([.5 3.5 -.3 .3])
  
  r = roots(coef)'
  plot(r,0*r,'ro')
  
  Jfd = full(au_jacobian_fd(@roots, coef, [], 1e-7));
  [r,J] = au_grad_roots(coef);
  au_prmat(J, Jfd, J./Jfd(:,2:end));
  return
end

if 0
  %% Derivation
  syms a b c d x real
  % Make polynomial with roots a,b,c,d and leading coeff 1
  p = collect(expand((x-a)*(x-b)*(x-c)*(x-d)), x)
  coef = fliplr(expand(coeffs(p, x)))
  au_assert_equal coef(1) 1  
  
end

if all(sort(size(coefs)) == [1 5])
  au_assert_equal coefs(1) 1
  r = roots(coefs);
  a = r(1);
  b = r(2);
  c = r(3);
  d = r(4);
  Jinv = [
[                -1,                -1,                -1,                -1]
[         b + c + d,         a + c + d,         a + b + d,         a + b + c]
[ - b*c - b*d - c*d, - a*c - a*d - c*d, - a*b - a*d - b*d, - a*b - a*c - b*c]
[             b*c*d,             a*c*d,             a*b*d,             a*b*c]
    ];
  svd(Jinv)
  J = inv(Jinv);
else
  error('unhandled');
end
