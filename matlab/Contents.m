%
% A collection of awf utilities
% 
% 
% The au_ prefix is because it's important in matlab's flat namespace
% that clashes of function names are avoided.
%
%
% Symbolic toolbox helpers
% 
% AU_COEFF  Coefficients of symbolic polynomial  
% AU_CCODE  Convert symbolic expression to C with common subexpression elimination  
% 
% 
% Faster/more convenient alternatives to matlab builtins
% 
% AU_SPARSE   Create spare matrices more efficiently  
% AU_WHIST    Weighted histogram  
% AU_BSX      Cleaner interface to bsxfun
% 
% 
% Printing and testing
% 
% AU_PRMAT   Compact matrix printing  
% AU_TEST*   Utilities for writing unit tests  
% AU_ASSERT* Easier assertions  
% 
% 
% Optimization/Numerics
% 
% AU_LEVMARQ     Simple Levenberg-Marquardt  
% AU_RANSAC      Simple RANSAC  
% AU_ROSENBROCK  Rosenbrock test function  
% AU_LOGSUMEXP   Efficient and safe log(sum(exp(V)))  
%
%
% MEX helper
% The library also contains au_mex.h, a clean C++ interface for
% mex files.   Type "edit au_mex.h" for more information.
