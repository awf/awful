function [e, J] = ellipse_err_lm(params, pts)
% This function is defined in terms of a function ellipse_ad_levmarq,
% as follows
% Call this function f, it is called with args thus:
%  f([p t1 ... tn]);    % p is 5 ellipse parameters, n is number of data
% Call ellipse_ad_mkpts g, it is called for each datum as
%  g(P) = g([p ti])
% So f is defined as
%  f(2*(i-1)+[1:2]) = g([p ti])
%

% extract params
geom = params(1:5);
tvalues = params(6:end);
n = length(tvalues);

% codegen from ellipse_ad_levmarq.m
% See  autogen_ellipse_ad_levmarq.cxx
dJ = autogen_ellipse_ad_levmarq([repmat(geom(:), [1 n]); tvalues(:)'], pts, true);
% au_assert('size(dJ,1)==14');

d = dJ(1,:);
e = d';

if nargout > 1
  Jgeom = dJ(2:6,:);
  Jt = dJ(7,:);
  J = sparse(2*n, 5+n);
  J(:,1:5) = Jgeom';
  cols = 5+(1:n);
  cols = [cols; cols];
  cols = cols(:);
  J = J + sparse(1:2*n, cols, Jt);
end
