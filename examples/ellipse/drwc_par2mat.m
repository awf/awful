function A = CONPAR2MAT(par)
% CONPAR2MAT      Convert conic from CONPAR to matrix (CONMAT) form.

% Author: Andrew Fitzgibbon, Edinburgh University AI Dept.
% Email: andrewfg@ed.ac.uk
% Date: 04 Apr 96

if min(size(par)) ~= 1
  error ('It is hard to store the conversion of several conics to matrix form at once.')
end
Ao = par(1);
Ax = par(2);
Ay = par(3);
Axx = par(4);
Ayy = par(5);
hAxy = par(6)/2;

A = [Axx hAxy Ax/2
    hAxy  Ayy Ay/2
     Ax/2 Ay/2 Ao];
