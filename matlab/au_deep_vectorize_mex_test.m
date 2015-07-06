

%% Test

au_test_begin

a = 2;
b.a = randn(2,3,4);
b.b = [1 1];
c = {randn(3,1), randn(1,1,2)};
b.c.d = c;
d = cell(2,2,3);
b1 = b;
b1.b = [2 3];
obj1 = {a,{b,c;c,b1}};

N = 40000;
tic
for k=1:N
  x1 = au_deep_vectorize_mex(obj1);
end
t1 = toc/N;
tic
N = 400;
for k=1:N
  x2 = au_deep_vectorize(obj1);
end
t2 = toc/N;

SPEEDUP = t2/t1;

au_test_equal x1 x2

%
obj1b = au_deep_unvectorize_mex(obj1, x1);
au_test_equal obj1 obj1b

obj2 = struct;
obj2(1).a = 1;
obj2(2).a = 1.2;
obj2(3).a = 2;
x1 = au_deep_vectorize(obj2);
[obj2a, n] = au_deep_unvectorize(obj2, x1);

% fprintf('x len = %d\n', length(x1));
au_test_equal n length(x1)
au_test_equal obj2a obj2

sx = sym('x', size(x1));
au_test_assert isa(sx,'sym')
obj2b = au_deep_unvectorize(obj2, sx);
au_test_assert isa(obj2b(1).a,'sym')

au_test_end

