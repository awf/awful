#include "au_mex_unchecked.h" // about 10% faster to disable bounds checks

// Declare mlx_function (C++ version of mexFunction)
void mlx_function(mlx_inputs& in, mlx_outputs& out)
{
   mlx_array<mlx_double> A(mlx_size {1,1}, in[0]); // Get input 0
   
   int size = int(A[0]);

   mlx_make_array<double> sum(size, size); // Make output array

   // Perform the operation
   for(mwSize j = 0; j < size; ++j) {
     double base = j+2;
     for(mwSize i = 0; i < size; ++i) {
       sum(i,j) = base;
       base += 1.0;
     }
   }
   out[0] = sum; // Assign to output
}
