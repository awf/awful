function au_mexall
d = cd;
[newdir,~,~] = fileparts(mfilename('fullpath'));

try
    cd(newdir)
    domex -largeArrayDims au_sparse.cxx
    domex au_whist.cxx
    domex au_deep_vectorize_mex.cxx
    domex au_deep_unvectorize_mex.cxx
catch e
    cd(d)
    rethrow(e)
end    
cd(d)

function domex(varargin)
disp(au_reduce(@(x,y) [x ' ' y], varargin, 'mex'))
mex(varargin{:})
