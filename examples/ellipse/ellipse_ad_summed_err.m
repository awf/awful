
function e = ellipse_ad_summed_err(params, data)

% NN  Distance from point to point on ellipse
% See DEMO_GANDER for use.
% If an ellipse is the parametric function
%    C(params, t) = [x(params; t), y(params; t))
% this computes
%    norm(data - C(params(1:5), params(6)).^2
% where
%  params: [cx cy rx ry theta t]
%  data: [x; y]

% Resulting mex has opcounts as follows:
%  sincos    4
%  *       104
%  +-       47
%  /         0
%  :=       94 [upper bound]

% That computes function, 6x1 gradient, 6x6 Hessian

%%
d = ellipse_ad_levmarq(params, data);
e = sum(d.^2);

if 0
  %% To make code
  au_autodiff_generate(@ellipse_ad_summed_err, rand(6,1), rand(2,1), ...
    'autogen_ellipse_ad_summed_err.cxx', 'HESSIAN=1')
end
