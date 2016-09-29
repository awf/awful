function Aout = demo_draw_conic(choice)

% DEMO_DRAW_CONIC A function
%               ...

% Author: Andrew Fitzgibbon <awf@robots.ox.ac.uk>
% Date: 09 Aug 97

if nargin == 0
  disp('Press a key to cycle through example conics.');
  for k=0:11
    draw_conic_demo(k);
    pause;
  end
  return
end

% A matrix
if choice == 0
  A = drwc_geom2mat([1 1 .9 .9 0]);
elseif choice == 1
  A = drwc_geom2mat([1 1 1.1 1.1 0]);
elseif choice == 2
  A = drwc_geom2mat([2, 1.5, 2, 1.5, 0]);
elseif choice == 9
  A = drwc_geom2mat([1 1 1 .07 30]);
elseif choice == 4
  A = drwc_geom2mat([1 1 3 0.5 40]);
elseif choice == 5
 A = drwc_geom2mat([1, 1, .2, 3, 60]);
 elseif choice == 6
 A = drwc_geom2mat([1, 1, 100,0.03, 30]);
elseif choice == 7
 A = drwc_par2mat([.7 -.4 -2 0 1 0]);
elseif choice == 8 
  A= [ -2 2 0; 2 2 0.15 ; 0 0.15 -1];
elseif choice == 10
 A = drwc_geom2mat([1 1.5 .05 -50 30]);
elseif choice == 12
 A = drwc_par2mat([2 -2 -6 1 5 2]);
elseif choice == 3
  A = drwc_geom2mat([1, 1, -1, 1, 60]);
elseif choice == 11
  A = drwc_geom2mat([1, 1, -1e-16 1e-16, 60]);
end

hold off
plot(nan)
hold on
axis([-1 3 -1 3]);

bbox = [0 2 0 2];
plot(bbox([1 2 2 1 1]), bbox([3 3 4 4 3]),'r');

draw_conic(A, bbox);

if nargout > 0
  Aout = A;
end
