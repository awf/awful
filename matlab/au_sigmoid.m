function s=  au_sigmoid(slope, x)

% AU_SIGMOID   s = 1./(1 + exp(-4*slope*x))
%               S = AU_SIGMOID(K, D)

if nargin == 0
  %% Demo
  t = -2:.01:2;
  hold off
  plot(t, au_sigmoid(1.0,t))
  hold on
  plot(t,0*t+1,'k');
  plot(t,0*t,'k');
  axis equal
  return
end

s = 1./(1 + exp(-4*slope.*x));
