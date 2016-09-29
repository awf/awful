function curves = draw_conic(a,original_handle)
% DRAW_CONIC    Draw a conic section.
%               DRAWCONIC(A) draws the conic section defined by
%               parameter vector A = [Ao Ax Ay Axx Ayy Axy]:
%               Axx x^2 + Axy x y + Ayy y^2 + Ax x + Ay y + Ao = 0

%               DRAWCONIC(H,HANDLE) alters the [XY]Data of HANDLE.
%               DRAWCONIC(H,'style') sets the linestyle.
%
%               Apologies for the peculiar ordering of the parameter
%               vector -- it makes some sense, as shorter vectors draw the
%               appropriate things:
%                 draw_conic([a b c])           line,
%                 draw_conic([a b c d])         circle,
%                 draw_conic([a b c d e])       unrotated conic. 

[x,y] = old_dc_pts(a);

%% Actually draw it
if nargin == 2,
  if isstr(original_handle)
    % it's a style
    curves = plot(x,y, original_handle);
  else
    set(original_handle, 'xdata', x, 'ydata', y);
    curves = original_handle;
  end
else
  curves = plot(x,y);
end
