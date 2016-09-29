function geom = drwc_par2geom(par)
% CONPAR2GEOM      Convert [Ao..Axy] conic to [cx cy rx ry theta].
%            Theta is in degrees.

% Author: Andrew Fitzgibbon, Edinburgh University AI Dept.
% Email: andrewfg@ed.ac.uk
% Date: 04 Apr 96

if nargin < 1
  %% test
  c = [.1 .3 .5 .7 31];
  p = drwc_geom2par(c);
  au_test_equal drwc_par2geom(p) c 1e-8
  return
end

thetarad = 0.5*atan2(par(:,6),par(:,4) - par(:,5));
cost = cos(thetarad);
sint = sin(thetarad);
sin_squared = sint.*sint;
cos_squared = cost.*cost;
cos_sin = sint .* cost;

Ao = par(:,1);
Au =   par(:,2) .* cost + par(:,3) .* sint;
Av = - par(:,2) .* sint + par(:,3) .* cost;
Auu = par(:,4) .* cos_squared + par(:,5) .* sin_squared + par(:,6) .* cos_sin;
Avv = par(:,4) .* sin_squared + par(:,5) .* cos_squared - par(:,6) .* cos_sin;

% ROTATED = [Ao Au Av Auu Avv]

tuCentre = - Au./(2.*Auu);
tvCentre = - Av./(2.*Avv);
wCentre = Ao - Auu.*tuCentre.*tuCentre - Avv.*tvCentre.*tvCentre;

uCentre = tuCentre .* cost - tvCentre .* sint;
vCentre = tuCentre .* sint + tvCentre .* cost;

Ru = -wCentre./Auu;
Rv = -wCentre./Avv;

Ru = sqrt(abs(Ru)).*sign(Ru);
Rv = sqrt(abs(Rv)).*sign(Rv);

geom = [uCentre, vCentre, Ru, Rv, thetarad * (180/pi)];
