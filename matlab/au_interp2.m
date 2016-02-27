%VGG_INTERP2  Fast 2d interpolation for images
%
%	V = vgg_interp2(A, X, Y)
%	V = vgg_interp2(A, X, Y, interp_mode)
%	V = vgg_interp2(A, X, Y, interp_mode, oobv)
%
% 2d interpolation on a regular grid - similar to matlab's interp2() but
% faster, and supports multiple channels and types. Note
% that while results of 'linear' and 'nearest' interpolation are the same
% as those of interp2(), those of cubic are not - vgg_interp2 uses a cubic
% hermite spline that is very fast to compute, unlike the natural cubic
% spline employed by interp2(), which does, however, yield a smoother
% interpolation.
%
%IN:
%	A - HxWxC double, single, uint16, uint8 or Logical array.
%	X - MxN horizontal offsets (1 being the centre of the first pixel).
%	Y - MxN vertical offsets (1 being the centre of the first pixel).
%	interp_mode - string, either 'cubic', 'linear' or 'nearest'. Default:
%	              'linear'.
%	oobv - 1x1 Out of bounds value. Default: NaN.
%
%OUT:
%	V - MxNxC interpolated values. Class is the same as that of oobv.

% Authors: Yoni Wexler, Tomas Werner, ojw, awf

error('Compile au_interp2.cxx');


%% Derivations for cubic
syms u v real
u1 = u * u;
v1 = v * v;
u2 = u1 * u;
v2 = v1 * v;
syms c0 c1 c2 c3 real
for m=0:3
  a = (c3 + c1) - (c2 + c0);
  bm = v2 * a + v1 * ((c0 - c1) - a) + v * (c2 - c0) + c1;
  dbmdv = 3*v1*a + 2*v*((c0 - c1) - a) + (c2 - c0);
end
a = (b(3) + b(1)) - (b(2) + b(0));
Bj = u2 * a + u1 * (b(0) - b(1) - a) + u * (b(2) - b(0)) + b(1);
