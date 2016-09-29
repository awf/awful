function [ex, ey] = draw_conic_evolute(a)

% DRAW_CONIC_EVOLUTE A function
% 

% Author: Andrew Fitzgibbon <awf@robots.ox.ac.uk>
% Date: 31 Oct 96

% Andrew W. Fitzgibbon, Artificial Intelligence, Edinburgh Univ.
% April 95

if size(a) == [3 3]
  A = a;
  if max(max(abs(a - a'))) > eps
    error('A not symmetric');
  end
  a = [A(3,3) A(1,3)*2 A(2,3)*2 A(1,1) A(2,2) A(1,2)*2];
end

% Copyright (C) Andrew Fitzgibbon, 1 Jan 1995.
% andrewfg@ed.ac.uk
% The University of Edinburgh
Ao = a(1);
Ax = a(2);
Ay = a(3);
Axx = a(4);
Ayy = a(5);
Axy = a(6);

% Is it a parabola?
discriminant = Axy^2 - 4*Axx*Ayy;
  
% Set up matrix form
A = [Axx Axy/2; Axy/2 Ayy]; b = [Ax Ay]';
  
% Generate normals
npts = 30;
theta = linspace(0,pi,npts);
u = [cos(theta);
  sin(theta)];
  
%%%% Solve for points at which gradient points along each normal
Ai = inv(A);
lambda = sqrt((b'*Ai*b - 4*Ao) ./ sum(u .* (Ai * u)));
  
% NaN out imaginary roots
I = find(imag(lambda)>eps);
lambda(I) = nan*I;
  
% Check for totally invisible conic
if all(isnan(lambda))
  % Disc is Axx Ax^2 + Ax*(Axy Ay) + Ayy Ay^2 - 4 C
  fprintf(2,'draw_conic2: no points (discriminant b''*A*b - 4 C = %g)\n', b'*A*b - 4*Ao)
end

% Continue to solve.
L = lambda([1 1],:); 
B = b(:,ones(1,npts));
Xpos = 0.5 * Ai * ( L.*u - B);
Xneg = 0.5 * Ai * (-L.*u - B);
X = [Xpos Xneg];
x = X(1,:);
y = X(2,:);

Gx = 2 * Axx * x +     Axy * y + Ax;
Gy =     Axy * x + 2 * Ayy * y + Ay;
G = [Gx Gy];

H = [
2*Axx    Axy;
  Axy  2*Ayy];

% Calc K

Gx2 = Gx.*Gx;
Gy2 = Gy.*Gy;

gradmag = sqrt(Gx2 + Gy2);

% Compute adjoint of H
J = inv(H) * det(H);

% grad(F)'*adj(H)*grad(F) = - \kappa |grad(F)|^3
K = -(Gx2 * J(1,1) + 2*Gx.*Gy * J(1,2) + Gy2 * J(2,2)) ./ (gradmag.^3);

R = 1 ./ K;

Nx = Gx ./ gradmag;
Ny = Gy ./ gradmag;

% plot(x + Nx.*R, y + Ny.*R, 'r');

%%% Faster:
Delta = -(Gx2 + Gy2) ./ (Gx2 * J(1,1) + 2*Gx.*Gy * J(1,2) + Gy2 * J(2,2));
plot(x + Gx .* Delta, y + Gy .* Delta, 'g');