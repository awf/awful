function e = ellipse_ad_levmarq(params, data)

% ELLIPSE_AD_LEVMARQ  Distance from point to point on ellipse
% See DEMO_GANDER for use.
% If an ellipse is the parametric function
%    C(params, t) = [x(params; t), y(params; t))
% this computes
%    data - C(params(1:5), params(6))
% where
%  params: [cx cy rx ry theta t]
%  data: [x; y]

% That computes 2x1 function, 12x1 gradient

%%
[xs,ys] = ellipse_points(params(1:5), params(6));
e = [xs; ys] - data(:);

if 0
  %% To make code
  au_autodiff_generate(@ellipse_ad_levmarq, rand(6,1), rand(2,1), 'autogen_ellipse_ad_levmarq.cxx')
end
