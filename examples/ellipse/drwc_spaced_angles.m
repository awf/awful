function angles = drwc_spaced_angles(phi1,phi2,l_infty)

% DRWC_SPACED_ANGLES Generate angles between phi1 and phi2,
%           weighted by distance to l_infty.

% Author: Andrew Fitzgibbon <awf@robots.ox.ac.uk>
% Date: 09 Aug 97

ang=[phi1];
i=1;
k=0.000005;
offset=0.1;
x=[cos(ang(i)),sin(ang(i)),1]';
dl_inf=l_infty'*x;
dl_inf2=dl_inf.*dl_inf;
gam=ang(i)+k*dl_inf2+offset;

while gam < phi2,
  i=i+1;
  ang=[ang,gam];
  x=[cos(ang(i)),sin(ang(i)),1]';
  dl_inf=l_infty'*x;
  dl_inf2=dl_inf.*dl_inf;
  gam=ang(i)+k*dl_inf2+offset;
end
angles=[ang,phi2];
