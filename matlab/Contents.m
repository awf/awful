function AUTO_GENERATED_FROM_Contents_template_txt
% The awf utility library: A collection of awf utilities
% The au_ prefix is because it's important in matlab's flat namespace
% that clashes of function names are avoided.
% 
% *Symbolic toolbox helpers*
% AU_COEFF        Extract polynomial coefficients from symbolic expr
% AU_CCODE        Generate optimized C code from symbolic expression.
% AU_AUTODIFF_GENERATE Generate code for function and derivatives
% 
% *Faster/more convenient alternatives to matlab builtins*
% AU_SPARSE       Create sparse matrices with low time/space overhead.
% AU_WHIST        Weighted histogram
% AU_BSX          A value class that implements a broadcastable data type
% 
% *Printing and testing*
% AU_PRMAT        Compact print of matrices.
% AU_DEEP_PRINT   Hierarchical print of object.
% AU_TEST*        Utilities for writing unit tests  
% AU_ASSERT*      Easier assertions  
% AU_ROSENBROCK   Rosenbrock
% AU_RUN_TESTS    Run all tests in the library
% au_progressbar_ascii A function
% 
% *Optimization*
% AU_OPTIMPROBLEM Simplified but powerful nonlinear least squares
% AU_LEVMARQ      Home-grown LM with line search
% AU_RANSAC       Ransac loop
% AU_DEEP_VECTORIZE Flatten arbitrary structure/cell a linear vector x.
% AU_DEEP_UNVECTORIZE Unflatten arbitrary structure/cell from a linear vector x.
% 
% *Special functions*
% AU_RODRIGUES    Convert axis/angle representation to rotation
% AU_LOGSUMEXP    Compute log(sum(exp(M))) stably
% AU_SIGMOID      s = 1./(1 + exp(-4*slope*x))
% 
% *File I/O helpers*
% AU_FSCAN_REGEXP File scan line by line splitting on regexp
% AU_STRIP_PATH   Remove directories matching REGEXP from PATH
% 
% *MEX helper*
% AU_MEX          C++ helper classes for MEX file writers.
% 
