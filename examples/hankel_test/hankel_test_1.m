function A = hankel_test_1(N)

% Generate a N-by-N matrix where A(i,j) = i + j;
A=zeros(N,N);
for jj = 1:N
  for ii = 1:N
    A(ii,jj) = ii + jj;
  end
end

