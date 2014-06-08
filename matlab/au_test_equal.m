function au_test_equal(expr1,expr2,tol,FORCE_PRINT)
% AU_TEST_EQUAL  Test all(EXPR1 == EXPR2), print result
%             au_test_equal det(M) 0 1e-7
%             au_test_equal('det(M)','0',1e-7);
%             au_test_equal 16+1 17
%             We call with strings rather than values to give much better
%             error messages.  The strings are evaluated in the caller's
%             context, so should behave sensibly.

if nargin == 0
  % Test
  au_test_test
  return
end

if nargin < 2
    error('au_test_equal: need at least two arguments')
end

if nargin < 3
  tol = 0;  % There is no other sensible default.
else
  if ischar(tol)
    tol = evalin('caller', tol);
    au_assert isnumeric(tol)
  end
end

if nargin < 4
  FORCE_PRINT = 0;
else
  if ischar(FORCE_PRINT)
    FORCE_PRINT = strcmp(FORCE_PRINT, 'print');
  end
end

mfile = au_mfilename(-1);
exprval1 = evalin('caller',expr1);
exprval2 = evalin('caller',expr2);
symbolic = isa(exprval1,'sym') || isa(exprval2,'sym');

hd = ['au_test_equal[' mfile ']:'];
err = inf;

if symbolic
    eq = all(exprval1 == exprval2);
elseif isequal(exprval1, exprval2)
    eq = 1;
elseif isnumeric(exprval1) && isnumeric(exprval2) && isequal(size(exprval1), size(exprval2))
    % Check for doubles within tolerance
    err = max(abs(double(exprval1(:)) - double(exprval2(:))));
    eq = err <= tol;
else
    eq = 0;
end

if ~eq
  if ~isnumeric(exprval1)
    fprintf(2, '%s\n', [hd ' FAILED: ' expr1 ' == ' expr2 ]);
  else
    fprintf(2, '%s *FAILED* |%s - %s| < %g (ERR=%g)\n', hd, expr1, expr2, tol, err);
    if (numel(exprval1) < 10 || FORCE_PRINT)
      m2s = @(x) mlp_mat2str(x, 3, 40);
      if ~strcmp(expr1, num2str(exprval1(:)'))
        fprintf(2, '   with %s = %s\n', expr1, m2s(exprval1));
      end
      if ~strcmp(expr2, num2str(exprval2(:)'))
        fprintf(2, '   with %s = %s\n', expr2, m2s(exprval2));
      end
    end
  end
else
  if tol ~= 0
    fprintf(1, '%s passed: |%s - %s| < %g\n', hd, expr1, expr2, tol);
  else
    fprintf(1, '%s passed: %s == %s\n', hd, expr1, expr2);
  end
end
