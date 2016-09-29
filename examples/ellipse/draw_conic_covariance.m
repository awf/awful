function draw_conic_covariance(x)

% DRAW_CONIC_COVARIANCE Draw conic corr to cov matrix of some pts
%               ...

% Author: Andrew Fitzgibbon <awf@robots.ox.ac.uk>
% Date: 02 Oct 03


mu = mean(x);
C = inv(cov(x));
Cu = C * mu';

A = [C -Cu; -Cu' (mu*Cu - .1)];

draw_conic(A);
