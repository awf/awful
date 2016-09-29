function [e, gradient, Hessian] = ellipse_err_ad(params, pts, conic_handle)
% This function is defined in terms of a function ellipse_ad_mkpts,
% as follows
% Call this function f, it is called with args thus:
%  f([p t1 ... tn]);    % p is 5 ellipse parameters, n is number of data
% Call ellipse_ad_mkpts g, it is called for each datum as
%  g(P) = g([p ti])
% So f is defined as
%  f = sum_i g([p ti])
% Now, we have access to g's 6x1 gradient and 6x6 Hessian,
%
%   dg/dP = [dg/dp]     d2g/dP2 = [ d2g/dp2   d2g/dtdp]
%           [dg/dt]               [ d2g/dpdt  d2g/dt2 ]
%
% So f's gradient (of size 5+n) is
%     the *sum* of the p gradients,
%     followed by a *list* of the t gradients
%  df/dp = sum_i dg{i}/dp
%  df/dti =  dg{i}/dt
%
% And the Hessian is [sum, list of sum; list of sum; list]
%  d2f/dp2 = sum_i d2g{i}/dp2
%  d2f/dpdti =  d2g{i}/dtdp
%  d2f/dti2 =  d2g{i}/dt2

% extract params
geom = params(1:5);
tvalues = params(6:end);
n = length(tvalues);

dJH = autogen_ellipse_ad_summed_err([repmat(geom(:), [1 n]); tvalues(:)'], pts, 2);
d = dJH(1,:);
J = dJH(2:7,:);
H = dJH(8:end,:);
e = sum(d);
gradient = zeros(5 + n, 1);
gradient(1:5) = sum(J(1:5,:), 2);
gradient(6:end) = J(6,:);

if nargout > 2
  H11 = zeros(5,5);
  H66 = zeros(6,6);
  HB = zeros(n,5);
  Hdiag = zeros(n,1);
  inds = au_tril_indices(6);
  for i=1:n
    % Hessian for [params, ti] is of form
    % [ A  b ]
    % [ b' c ]
    %
    H66(inds) = H(:,i);
    H66s = (H66 + H66')./(1 + eye(6));

    H11 = H11 + H66s(1:5,1:5);
    b = H66s(1:5,6);
    c = H66s(6,6);
    HB(i,:) = b';
    Hdiag(i) = c;
  end
  ii = 1:n;
  Hessian = sparse([],[],[],5+n,5+n,11*n+25);
  Hessian(1:5, 1:5) = H11;
  Hessian(5+ii, 1:5) = HB;
  Hessian(1:5, 5+ii) = HB';
  Hessian = Hessian + sparse(5+ii,5+ii, Hdiag);
end

if nargin > 2
  draw_conic(g2p(geom), conic_handle);
  drawnow
end
end
