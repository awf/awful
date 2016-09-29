function [x,y,t] = au_conic_closest_point(params, points, fast)

% CONIC_CP2     Conic closest point.
%               PARAMS = [Ao Ax Ay Axx Ayy Axy]
%               [x,y] = au_conic_closest_point(params, points[, FAST])
%               distances_squared = au_conic_closest_point(params, points)
%               FAST=1 disablses asserts
%               

%% TEST 
if nargin == 0
  %%
  pp=[    1.7515   -1.4221   -2.0054    0.3620    0.6766    0.6877 ];
  
  pts = [1 1]';
  au_conic_closest_point(pp, pts)
  return
  
  params = drwc_geom2par([.6 .5 .4 .2 23]);
  hold off
  plot(nan)
  axis([0 1 0 1]);
  hold on
  draw_conic(params);
  
  pts = rand(2,10);
  plot(pts(1,:), pts(2,:), '.');
  [cpx,cpy] = au_conic_closest_point(params, pts);
  plot(cpx, cpy, 'r.')
  plot([pts(1,:)' cpx']', [pts(2,:)' cpy']', 'r-')
  axis equal
  
  
end

if 0
  %% Derivation
  syms a b c s u v lambda real
  % Find t to minimize norm([a cos t, b sin t] - [u v])
  % Do it by finding [c s] on unit circle.
  e = (a*c - u)^2 + (b*s - v)^2 - lambda*(s^2 + c^2 - 1)
  G = gradient(e, [s c lambda])
  lambda_from_1 = solve(G(1), lambda)
  SG2 = simplify(subs(G(2), lambda, lambda_from_1))
  % 2*b*c*(y - b*s) - 2*a*s*(x - a*c)
  EqC = subs(G(3), s, solve(SG2, s))
  EqC = collect(EqC * ((- c*a^2 + u*a + c*b^2)^2), c)
  coef = au_coeff(EqC, c)'
  
end

%% Params
if nargin < 3
  fast= 0;
end

%% Code
Ao  = params(1);
Ax  = params(2);
Ay  = params(3);
Axx = params(4);
Ayy = params(5);
Axy = params(6);

P = points;
if ~fast, au_assert_equal('size(P,1)', '2'); end

% Rotate the conic
theta = 0.5*atan2(Axy,Axx - Ayy);
cost = cos(theta);
sint = sin(theta);
sin2 = (sint)^2;
cos2 = (cost)^2;

% rotation-free
Auu = Axx * cos2 + Ayy * sin2 + Axy * sint * cost;
Avv = Axx * sin2 + Ayy * cos2 - Axy * sint * cost;
Au =   Ax * cost + Ay * sint;
Av = - Ax * sint + Ay * cost;

% calculate translation
tu = Au/(2*Auu);
tv = Av/(2*Avv);
C = Ao - Auu*tu*tu - Avv*tv*tv;
Auu = -Auu / C;
Avv = -Avv / C;

X0 = P(1,:);
Y0 = P(2,:);

% transform start points
U0 =  cost * X0 + sint * Y0 + tu;
V0 = -sint * X0 + cost * Y0 + tv;

% temps
U02 = U0.*U0;   V02 = V0.*V0;
ai = 1 / Auu;  bi = 1 / Avv;

% polynomial coeffs
N = size(P,2);
C = ones(5,N);
l = ones(1,N);
C(1,:) = l;
C(2,:) = -2.0*(ai+bi)*l;
C(3,:) = -(U02*ai+V02*bi-ai*ai-4.0*ai*bi-bi*bi);
C(4,:) = 2.0*(U02+V02-ai-bi)*ai*bi;
C(5,:) = -(U02*bi+V02*ai-ai*bi)*ai*bi;

R = zeros(4,N);
Rimag = false(4,N);
for n = 1:N
  % Extract roots
  r = roots(C(:,n));
  R(:,n) = r;
  Rimag(:,n) = abs(imag(r)) > eps;
end

%% Calculate distances
U0pad = U0([1 1 1 1], :);
V0pad = V0([1 1 1 1], :);

Upad = U0pad ./ (1 - R*Auu);
Vpad = V0pad ./ (1 - R*Avv);
  
du = Upad - U0pad;
dv = Vpad - V0pad;
dists2 = du.*du + dv.*dv;

dists2(Rimag) = 1e20;

% check dists are finite
if ~fast
  au_assert('all(isfinite(dists2(:)))')
end

[dmin, imin] = min(dists2);

if nargout == 1
  x = dmin;
else
  %% Rotate point and gradient back onto conic
  inds = sub2ind(size(Upad), imin, 1:N);
  u = Upad(inds);
  v = Vpad(inds);
  u = u - tu; v = v - tv;
  x = cost * u - sint * v;
  y = sint * u + cost * v;
  
  if nargout == 3
    t = atan2(v,u);
  end
end
