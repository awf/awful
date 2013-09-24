function [x, unpack] = au_fmin_pack(varargin)
% AU_FMIN_PACK  Pack variously-sized parameters into a linear vector x.
%         [x, unpack] = au_fmin_pack(aVec,aMat,aScalar)
%         aVec == unpack(x,1)
%         aMat == unpack(x,2)
%         [aVec,aMat,aScalar] = unpack(x)

% awf, aug13

if nargin == 0
  % Test
  a = 2;
  b = randn(2,3,4);
  c = randn(3,1);
  [tx, tunpack] = au_fmin_pack(a,b,c);
  
  [a1,b1,c1] = tunpack(tx);

  au_test_equal tunpack(tx)   a
  au_test_equal tunpack(tx,2) b

  au_test_equal a a1
  au_test_equal b b1
  au_test_equal c c1

  return
end


sz = 0;
for k=1:nargin
  sz = sz + numel(varargin{k});
end

unpack_data(sz).start = [];
unpack_data(sz).sz = [];

x = zeros(sz,1);
start = 1;
for k=1:nargin
  Ak = varargin{k};
  szk = numel(Ak);
  x(start:start+szk-1) = Ak(:);
  unpack_data(k).start = start;
  unpack_data(k).sz = size(Ak);
  
  start = start + numel(varargin{k});
end

unpack = @(x, varargin) unpack_from_data(x, unpack_data, varargin{:});

function [varargout] = unpack_from_data(x, ud, k)
if nargin == 2
  todo = 1:nargout;
else
  todo = k;
end
varargout = cell(nargout,1);
for kk=1:length(todo)
  k=todo(kk);
  start = ud(k).start;
  sz = ud(k).sz;
  Ak = reshape(x(start:start+prod(sz)-1), sz);
  varargout{kk} = Ak;
end
