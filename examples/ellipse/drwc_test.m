function drwc_test(l_infty, d)

% DRWC_TEST     A function
%               ...

% Author: Andrew Fitzgibbon <awf@robots.ox.ac.uk>
% Date: 05 Sep 97

if nargin == 2
  % angle, distance
  theta = l_infty;
  l_infty = [cos(theta) sin(theta) d]; 
end

a = l_infty(1);
b = l_infty(2);
c = l_infty(3);

t = linspace(-pi,pi,64);

% Image of unit circle under homography [I_2x2 0_2x1; l_infty]:
% C(t) is homg parametrized curve
x = cos(t);
y = sin(t);
z = a * cos(t) + b * sin(t) + c;
C = [x;y;z];

% Derivatives
dx = -sin(t);
dy = cos(t);
dz = -a * sin(t) + b * cos(t);
dC = [dx;dy;dz];

% nc = nonhomogeneous c
ncx = x ./ z; 
ncy = y ./ z;

dnx = (x .* dz - z .* dx) ./ (z.*z);
dny = (y .* dz - z .* dy) ./ (z.*z);

%%% Start graphics

% Init axes
hold off
plot(nan)
hold on
axis([-2 2 -1 1] * 8)
axis equal

% Draw tangent lines
for k = 1:length(t)
  % tangent_line = cross(C(:,k), dC(:,k));
  tt = t(k);
  tangent_line = [a + c * cos(tt), b + c * sin(tt), -1];
  
  h = draw_implicitline(tangent_line); set(h, 'color', [1 1 1]/3);
end

% Plot sample points
% plot(ncx, ncy, 'yo')

% Draw conic
plot(ncx, ncy)

handles = [];
choose_l = 0;
while 1
  [lx, ly] = awfgetline; 
  L = cross([lx(1) ly(1) 1], [lx(2) ly(2) 1]);
  if norm(L) == 0
    disp('drwc_test: returning');
    return
  end

  L = L / norm(L(1:2));
  
  delete(handles)
  handles = [];
  handles = [handles draw_implicitline(L)];

  l1 = L(1);
  l2 = L(2);
  l3 = L(3);
  
  % cos_tmax = l1 * (c^2 + 1) / (c * (1 - 2*l1));
  % cos_tmax = ((l1 - 1) + sqrt(l1*(c*c + l1 - 1)))/c;
  lambda = l1 / l2;

  for s = [-1 1]
    if s == -1, col = 'g'; else, col = 'm'; end;

    % Choose cos_tmax:
    disc = c^2 - l2^2;
    if disc < 0
      disp('Imag soln')
    else
      c_cosphi = s * sqrt(disc);
      
      c_cos_tmax = l1 * c_cosphi - l2 * l2;
      c_sin_tmax = l2 * c_cosphi + l1 * l2;
      
      xmax = c_cos_tmax / c;
      ymax = c_sin_tmax / c;
      zmax = a * xmax + b * ymax + c;
      if zmax == 0
	disp('Soln at infinity');
      else
	p1 = [xmax, ymax]/zmax;
	
	% furthest point
	h2 = plot(p1(1),p1(2), 'ro');
      
	% tangent
	h3 = draw_implicitline([a + c * xmax, b + c * ymax, -1], 'b');
    
	D = s * l3 + 1 / sqrt(c^2 + 2*c_cos_tmax + 1);
    
	
	n = orth([a + c * xmax; b + c * ymax]);
	h4 = draw_line(p1, p1 - n' * D, col);
	handles = [handles h2 h3 h4];
      end
    end
  end
end
