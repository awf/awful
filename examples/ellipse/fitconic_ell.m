function [pvec,N,M] = fitconic_ell(P, cd, realpvec)
% Create design matrix

x=P(:,1);
y=P(:,2);

n = length(x);
O = x*0;
l = O+1;

D  = [ x.*x x.*y y.*y x y l ];

%
% Create scatter matrix
%
S = D'*D;

%
% Create constraint matrix: (2 Axx Ayy - Axy^2 + 2 Ayy Axx = 1)
%
C = [ 0  0  2  0  0  0
      0 -1  0  0  0  0
      2  0  0  0  0  0
      0  0  0  0  0  0
      0  0  0  0  0  0
      0  0  0  0  0  0 ];

pvec = cmin(S, -C);

% Swap components into [1 x y xx yy xy] order
pvec = [ pvec(6) pvec(4) pvec(5) pvec(1) pvec(3) pvec(2)];
if pvec(1) < 0, pvec = -pvec; end
