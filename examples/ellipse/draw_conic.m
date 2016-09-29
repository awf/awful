function curves = draw_conic(A,original_handle)
% DRAW_CONIC    Draw a conic section.
%               DRAWCONIC(A) draws the conic section defined by
%               parameter vector A = [Ao Ax Ay Axx Ayy Axy]:
%               Axx x^2 + Axy x y + Ayy y^2 + Ax x + Ay y + Ao = 0
%
%               DRAWCONIC(H,HANDLE) alters the [XY]Data of HANDLE.
%               DRAWCONIC(H,'style') sets the linestyle.
%
%               Apologies for the peculiar ordering of the parameter
%               vector -- it makes some sense, as shorter vectors draw the
%               appropriate things:
%                 draw_conic([a b c])           line,
%                 draw_conic([a b c d])         circle,
%                 draw_conic([a b c d e])       unrotated conic. 

% Authors: Manish Jethwa and Andrew Fitzgibbon, Robotics Research, Univ Oxford.
% Email: u94mj@robots.ox.ac.uk, awf@robots.ox.ac.uk
% Date: July 1997

% Convert to matrix form
if ~all(size(A) == [3 3])
  if min(size(A)) ~= 1
    error('draw_conic: A must be 3x3 or 1x[3-6]');
  end
  a = [A(:)' 0 0 0];
  a = a(1:6);
  A = [
    a(4)   a(6)/2 a(2)/2;
    a(6)/2 a(5)   a(3)/2
    a(2)/2 a(3)/2 a(1)
    ];
end

%% Actually draw it
if nargin == 2,
  if isstr(original_handle)
    % it's a style, clip to axis
    [x,y] = draw_conic_pts(A, axis);
    curves = plot(x,y, original_handle);
  elseif sort(size(original_handle)) == [1 4]
    % bbox supplied
    [x,y] = draw_conic_pts(A, original_handle);
    curves = plot(x,y);
    
  else
    % handle supplied, use bounds of parent axis
    ax = get(original_handle,'parent');
    bbox = [get(ax,'xlim') get(ax,'ylim')];
    [x,y] = draw_conic_pts(A, bbox);
    set(original_handle, 'xdata', x, 'ydata', y);
    curves = original_handle;
  end
else
  [x,y] = draw_conic_pts(A, axis);
  curves = plot(x,y);
end
