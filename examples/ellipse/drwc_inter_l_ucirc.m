function [theta1, theta2] = drwc_inter_l_ucirc(l)
% DRWC_INTER_L_UCIRC Intersect a line and the unit circle, homogeneous.
%               [theta1, theta2] = inter_line_circle(l);

a2=l(1)*l(1);
b2=l(2)*l(2);
c2=l(3)*l(3);
if c2<=(a2+b2),
  alpha=acos(sqrt(c2/(a2+b2)));
  if l(3) < 0,
    l = -l; 
  end
  beta=atan2(l(2),l(1))+pi;
  theta1=beta-alpha;
  theta2=beta+alpha;
else
  theta1=i;
  theta2=0; % fixme
end
