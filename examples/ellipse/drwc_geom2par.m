function par = drwc_geom2par(geom)
% CONGEOM2PAR      Convert conic from geometric to parameter-vector form.

% Author: Andrew Fitzgibbon, Edinburgh University AI Dept.
% Email: andrewfg@ed.ac.uk
% Date: 04 Apr 96

% in is [cx cy rx ry theta]
cx = geom(1);
cy = geom(2);
rx = geom(3);
ry = geom(4);
theta = geom(5)*pi/180;

Rx=sign(rx)/(rx^2);
Ry=sign(ry)/(ry^2);

theta = -theta;
cost = cos(theta);
sint = sin(theta);
sin2t = 2*sint*cost;
Axx = cost^2*Rx + sint^2*Ry;
Axy = sin2t*(Ry-Rx);
Ayy = Rx*sint^2 + Ry*cost^2;
Ax = sin2t*(Rx - Ry)*cy - 2*cx*(Ry*sint*sint + Rx*cost*cost);
Ay = sin2t*(Rx - Ry)*cx - 2*cy*(Ry*cost*cost + Rx*sint*sint);
Ao = Ry*(cx*sint + cy*cost)^2 + Rx*(cx*cost - cy*sint)^2 - 1;

par = [Ao Ax Ay Axx Ayy Axy];
