function au_assert(expr)
% AU_ASSERT  Assert all(EXPR), print expr if not
%             au_assert('det(M) > 0');

exprval = evalin('caller',expr);
if ~all(exprval(:))
  error(['au_assert: FAILED: ' expr]);
end
