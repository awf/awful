function [I, Scale] = au_colorize_signed_image(v, colorscheme, gamma)

% AU_COLORIZE_SIGNED_IMAGE Make an RGB image of a signed matrix.
%               Red is +ive, Blue -ive, Black/White for zero
%               I = au_colorize_signed_image(M);
%
%               I = au_colorize_signed_image(M, COLORSCHEME)
%               colorscheme is one of:
%                  1, 'whitebg' : red-white-blue handy for printouts
%                  2, 'blackbg' : red-black-blue for onscreen display
%                  3, 'gray'    : Classic grayscale, zero is mid-gray
%
%               I = au_colorize_signed_image(M, COLORSCHEME, GAMMA)
%               Apply gamma-correction factor to the image, useful to
%               accentuate small values, e.g. GAMMA = 0.2;
%
%               [I, Scale] = au_colorize_signed_image(...) returns the
%               scale factor max(abs(M(:)) which was applied to the data.
%               
%               In comparison to POLARMAP, this function actually returns
%               a new 8-bpp RGB image with the transformed values.  The
%               advantage is that you can save the image to a PNG for use
%               in documents, the disadvantage is that COLORBAR does not
%               do anything.

% Author: Andrew Fitzgibbon <awf@fitzgibbon.ie>
% Date: 07 Aug 2012

if nargin == 0
  F = peaks(256);
  F = sign(F).*F.^2;
  image([...
    au_colorize_signed_image(F) ...
    au_colorize_signed_image(F,'blackbg') ...
    au_colorize_signed_image(F,'gray');
    au_colorize_signed_image(F, 'whitebg', 0.5) ...
    au_colorize_signed_image(F,'blackbg', 0.25) ...
    au_colorize_signed_image(F,'gray', 0.25) ...
    ])
  axis image
  title('Demo of awf\_colorize\_signed\_image');
  return
end

if nargin < 2
  colorscheme = 1;
end

if nargin < 3
  gamma = 0;
end

vmax = max(v(:));
vmin = min(v(:));
vabsmax = max(abs(v(:)));
% fprintf('vmax = %g\n', vabsmax);

if vmax == 0
  intpos = zeros(size(v));
else
  intpos = (v > 0) .* v / vabsmax;
end

if vmin == 0
  intneg = zeros(size(v));
else
  intneg = -(v <= 0) .* v / vabsmax;
end

if gamma
  intneg = intneg.^gamma;
  intpos = intpos.^gamma;
end

switch colorscheme
case {1, 'whitebg'}
  % whitebg, handy for printouts
  I = uint8(255*cat(3, 1-intneg, v*0 + 1-intneg-intpos, 1 - intpos));
  
case {2, 'blackbg'}
  % blackbg
  I = uint8(255*cat(3, intpos, v*0, intneg));
case {3, 'gray'}
  % Classic grayscale
  I = repmat(uint8(v/vabsmax * 128 + 128), [1, 1, 3]);
otherwise
  error('colorscheme')
end

if nargout > 1
  Scale = vabsmax;
end
