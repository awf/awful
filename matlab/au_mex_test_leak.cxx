
// Bare-bones example of au_mex

#include "au_mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   mlx_make_array<double> tmp(1000,1000); // Make possible leaker
}
