#include "au_mex.h"

// Declare mlx_function (C++ version of mexFunction)
void mlx_function(mlx_inputs& in, mlx_outputs& out)
{
   mlx_array<double> A(in[0]); // Get input 0
   mlx_array<double> B(in[1]); // Get input 1

   mlx_make_array<double> sum(A.size); // Make output array

   // Perform the operation
   for(int i = 0; i < A.numel(); ++i)
     sum[i] = A[i] + B[i];

   out[0] = sum; // Assign to output
}
