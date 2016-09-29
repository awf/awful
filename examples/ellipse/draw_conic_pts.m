function [out_x,out_y] = draw_conic_pts(A, bbox)

% DRAW_CONIC_PTS Generate points for draw_conic (q.v.)
%              [x,y] = draw_conic_pts(A, bbox);

% Authors: Manish Jethwa and Andrew Fitzgibbon, Robotics Research, Univ Oxford.
% Email: u94mj@robots.ox.ac.uk, awf@robots.ox.ac.uk
% Date: July 1997

if max(max(abs(A - A'))) > 10*eps
  error('A not symmetric');
end

% Get bounding box
a_xmin = bbox(1); 
a_xmax = bbox(2); 
a_ymin = bbox(3); 
a_ymax = bbox(4);

% Check for rank-degeneracy
[R,D] = drwc_eigsrt(A);
rank_deficiency = sum(D == 0);
if rank_deficiency > 0
  if ~any(A(1:2,1:2))
    Ax = A(1,3)*2;
    Ay = A(2,3)*2;
    Ao = A(3,3);
    % Draw a line
    if abs(Ax) > abs(Ay),
      out_y = [a_ymin,a_ymax];
      out_x = (-1/Ax)*(Ay*out_y+Ao);
    else
      out_x = [a_xmin,a_xmax];
      out_y = (-1/Ay)*(Ax*out_x+Ao);
    end
  else
    disp('draw_conic: Denegerate conic: D =');
    disp(D)
  end
  return
end

% Define homography matrix H and inverse !
qd = sqrt(abs(D));
Q = diag(qd);
Qinv = diag(1./qd);

H=Q*R';
Hinv = R * Qinv;

% create an array containing the clipping lines
l = [1 0 -a_xmin]';
r = [1 0 -a_xmax]';
b = [0 1 -a_ymin]';
t = [0 1 -a_ymax]';

boundary_lines=[l,b,r,t];

% intersect preimages of lines with unit circle and find acceptable alpha values
[alphas,l_infty] = drwc_find_angles(boundary_lines,Hinv);

out_x = [];
out_y = [];

% For each segment
for m=1:2:length(alphas),

  % Generate samples such that output is evenly spaced.
  betas = drwc_spaced_angles(alphas(m),alphas(m+1),l_infty);
  circle_points = [
    cos(betas);
    sin(betas);
    ones(1,length(betas))
    ];

  % Transform circle points into original frame
  Q = Hinv * circle_points;
  
  % Add to output
  out_x = [out_x Q(1,:)./Q(3,:) nan];
  out_y = [out_y Q(2,:)./Q(3,:) nan];
end
