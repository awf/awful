function x = rand_donut(Rx,Ry,n)
% RAND_DONUT Return a set of points drawn from an ellipse + Gaussian
%             x = rand_donut(Rx,Ry,n);
%            To have a non-axis-aligned ellipse, rotate and translate x

t = rand(n,1) * 2*pi;
x = [Rx*cos(t) Ry*sin(t)];
x = x + randn(n,2);
