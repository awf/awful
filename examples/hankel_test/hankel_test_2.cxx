
// It's about 10% faster on this test to disable bounds checks
#include "au_mex_unchecked.h"

// Declare mlx_function (C++ version of mexFunction)
void mlx_function(mlx_inputs& in, mlx_outputs& out)
{
   mlx_array<mlx_double> A(mlx_size {1,1}, in[0]); // Get input 0
   
   int size = int(A[0]);

   mlx_make_array<double> sum(size, size); // Make output array

   // Perform the operation
   for(mwSize j = 0; j < size; ++j)
     for(mwSize i = 0; i < size; ++i)
       sum(i,j) = double(i+j+2);

   out[0] = sum; // Assign to output
}
