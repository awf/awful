function krondemo_test(A,B,H)

if nargin == 0
  evalin('base', 'krondemo')
  return
end

ntrials = 10;

tic
H = zeros(size(H));
for k=1:ntrials
  P = kron(A,B);
  H = H + P;
end
toc

tic
H2 = zeros(size(H));
for k=1:ntrials
  for j = 1:size(A,2)
    jrng = (j-1)*size(B,2)+[1:size(B,2)];
    for i = 1:size(A,1)
      irng = (i-1)*size(B,1)+[1:size(B,1)];
      P = A(i,j)*B;
      X = H2(irng,jrng);
      X = X  + P;
      H2(irng,jrng) = X;
    end
  end
end
toc

tic
H2a = zeros(size(H));
for j = 1:size(A,2)
  jrng = (j-1)*size(B,2)+[1:size(B,2)];
  for i = 1:size(A,1)
    irng = (i-1)*size(B,1)+[1:size(B,1)];
    X = H2a(irng,jrng);
    for k=1:ntrials
      P = A(i,j)*B;
      X = X  + P;
    end
    H2a(irng,jrng) = X;
  end
end
toc

tic
H3 = zeros(size(H));
for k=1:ntrials
  for j = 1:size(A,2)
    jrng = (j-1)*size(B,2)+[1:size(B,2)];
    Htmp = H3(:,jrng);
    for i = 1:size(A,1)
      irng = (i-1)*size(B,1)+[1:size(B,1)];
      P = A(i,j)*B;
      X = Htmp(irng,:);
      X = X  + P;
      Htmp(irng,:) = X;
    end
    H3(:,jrng) = Htmp;
  end
end
toc

au_assert_equal H H3
