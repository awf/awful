
// au_mex_example_2
// Shows how to specialize for different types.

#include "au_mex.h"

template <class Real>
mlx_cast<Real> Compute(mlx_cast<Real> const& A, mlx_cast<Real> const& B)
{
   mlx_make_array<Real> sum(A.size); // Make output array

   // Perform the operation
   for(int i = 0; i < A.numel(); ++i)
     sum[i] = A[i] + B[i];

   return sum;
}

template <class Real>
bool try_cast(mxArray const* pA, mxArray const* pB, mlx_output* out)
{
   mlx_cast<Real> A(pA);
   mlx_cast<Real> B(pB);
   if (!(A && B))
       return false;
   
   mlx_assert(A.size == B.size);// Check sizes match
   
   *out = Compute(A, B);
   return true;
}   

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   mlx_inputs  in(nrhs, prhs); // Wrap inputs
   mlx_outputs out(nlhs, plhs); // Wrap outputs
   
   // Enumerate the types.  You really do have to do this, so that the 
   // C++ compiler can lay down different code for each case.
   // You could clean this up with a macro if you like that sort of thing.
   if (try_cast<mlx_double>(in[0], in[1], &out[0])) return;
   if (try_cast<mlx_single>(in[0], in[1], &out[0])) return;
   if (try_cast<mlx_uint8>(in[0], in[1], &out[0])) return;
 
   mexErrMsgIdAndTxt("awful:types", "We don't handle this input type combo: %s, %s", 
           mxGetClassName(in[0]), mxGetClassName(in[1]));
}
