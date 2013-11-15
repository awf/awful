// Bare-bones example of au_mex

#include "au_mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   mlx_inputs  in(nrhs, prhs); // Wrap inputs
   mlx_outputs out(nlhs, plhs); // Wrap outputs

   mlx_cast<double> A(in[0]); // Get input 0
   mlx_cast<double> B(in[1]); // Get input 1

   mlx_assert(A.size == B.size);// Check sizes match

   mlx_make_array<double> sum(A.size); // Make output array

   // Perform the operation
   for(int i = 0; i < A.numel(); ++i)
     sum[i] = A[i] + B[i];

   out[0] = sum; // Assign to output
}
