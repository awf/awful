%% Demonstrate graph-coloring for finite difference Jacobian calculation
%

a = ones(3,2);
J = [blkdiag(a,a,a,a) ones(12,3)];

n = size(JacobPattern, 2);
p = colamd(JacobPattern)';
p = (n+1)*ones(n,1)-p;
group = color(JacobPattern,p);

hold off
colors = [
  1 0 0
  0 0 1
  0 0 0
  0 .7 0
  .7 .7 0
  ];
colors = {
  'r'
  'k'
  'r+'
  'b+'
  'b'
  };

spy(J)
 
%%
for k=1:max(group)
  spy(J.*au_bsx(group'==k), colors{k})
  hold on
end

%%
for k=1:d
  J(:,k) = F(x + e_k);
end
