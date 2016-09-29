function [x,y] = ellipse_points(params, t)

% sample points on ellipse.

cx = params(1);
cy = params(2);
rx = params(3);
ry = params(4);
theta = params(5);

x0 = rx * cos(t);
y0 = ry * sin(t);
c = cos(theta);
s = sin(theta);
x = cx + c * x0 - s * y0;
y = cy + s * x0 + c * y0;
