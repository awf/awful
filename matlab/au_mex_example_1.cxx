//#define AU_MEX_UNCHECKED /* define this to go fast */
#include "au_mex.h"

// Declare mlx_function (C++ version of mexFunction)
// Compare to 
// https://awful.codeplex.com/SourceControl/latest#matlab/au_mex_example_1.cxx 
void mlx_function(mlx_inputs& in, mlx_outputs& out)
{
   mlx_array<mlx_double> A(mlx_size{2,3}, in[0]); // Get input 0
   mlx_array<mlx_double> B(A.size, in[1]); // Get input 1

   mlx_make_array<double> sum(A.size); // Make output array

   // Perform the operation
   for(mwSize i = 0; i < A.numel(); ++i)
     sum[i] = A[i] + B[i];

   out[0] = sum; // Assign to output
}
