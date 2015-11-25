
mex -Ic:\dev\codeplex\awful\matlab -O hankel_test_2.cxx
mex -Ic:\dev\codeplex\awful\matlab -O hankel_test_3.cxx

report = @(msg,N) fprintf('%s: %.1fmsec\n', msg, toc/N*1000);

%%
N=2000
% Generate a N-by-N matrix where A(i,j) = i + j;
tic
for k=1:100
  A=zeros(N,N);
  for jj = 1:N
    for ii = 1:N
      A(ii,jj) = ii + jj;
    end
  end
end
report('Matlab JIT', k);

%%
tic
for k=1:100
  A2 = hankel_test_1(N);
end
report('Matlab JIT, function call', k);
au_assert_equal A A2

%%
broadcast = @au_bsx;
tic
for k=1:100
  A3 = broadcast(1:N) + (1:N)';
end
report('au_bsx', k);
au_assert_equal A A3

%%
tic
for k=1:100
  A4 = hankel_test_2(N);
end
report('mex -- straightforward', k);
au_assert_equal A A4

%%
tic
for k=1:100
  A5 = hankel_test_3(N);
end
report('mex -- double incremented', k);
au_assert_equal A A5
