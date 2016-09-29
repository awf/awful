function  [alphas,l_infty] = drwc_find_angles(b,Hinv)

% DRWC_FIND_ANGLES Intersect lines with unit circle.
%               ...

% Clipping boundary mappped onto unit circle
boundary_lines = Hinv'*b;

% Find preimage of line at infinity
l_infty = Hinv' * [0 0 1]';

% Collect angles of intersections of each clipping line with unit circle.
phi=[];
for i=1:4,
  [phi1,phi2]= drwc_inter_l_ucirc(boundary_lines(:,i));
  if phi1~=sqrt(-1), 
     while phi1<0,
       phi1 =phi1+2*pi;
     end
     while phi1>2*pi,
       phi1 =phi1-2*pi;
     end
     while phi2<0
       phi2 =phi2+ 2*pi;
     end 
     while phi2> 2*pi,
       phi2 = phi2 - 2*pi;
     end
     phi=[phi,phi1,phi2]; % place phi values in array
  end
end

% sort the accepted angles,
phi=sort(phi)+2*pi;
if isempty(phi)
  % Fully inside clipping region -- set start/end point to be
  % closest point to origin of l_infty
  % phi=[atan2(l_infty(2) / l_infty(3), l_infty(1) / l_infty(3))];
  if l_infty(3) < 0,
    s = -1;
  else
    s = 1;
  end
  phi=[atan2(l_infty(2) * s, l_infty(1) * s)];
end

% Wrap angles
phi=[phi,phi(1)+2*pi];

phis=[];
for j=1:length(phi)-1,
  alpha=(phi(j)+phi(j+1))/2;H_mp=Hinv*[cos(alpha),sin(alpha),1]';
  Hmx=H_mp(1)/H_mp(3);
  Hmy=H_mp(2)/H_mp(3);
  if (Hmx>-b(3,1))&(Hmx<-b(3,3))&(Hmy>-b(3,2))&(Hmy<-b(3,4)) % then accept
    phis=[phis,phi(j),phi(j+1)]; 
  end
end
alphas=phis;

