// B = vgg_interp2(A, X, Y[, method[, oobv]])

// $Id: vgg_interp2.cxx,v 1.4 2008/03/30 19:25:28 ojw Exp $
// Rewritten by ojw 20/9/06
// awf added derivatives, 12/02/16

#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "au_mex.h"

enum method_t {
  nearest,
  linear,
  cubic
};

template<class T, class U> 
void interp_nearest(mlx_array<U>& B, mlx_array<T>& A, mlx_array<mlx_double>& X, mlx_array<mlx_double>& Y, U oobv);

template<class T, class U> 
void interp_linear(mlx_array<U>& B, mlx_array<T>& A, mlx_array<mlx_double>& X, mlx_array<mlx_double>& Y, U oobv);

template<class T, class U> 
void interp_cubic(mlx_array<U>& B, mlx_array<T>& A, mlx_array<mlx_double>& X, mlx_array<mlx_double>& Y, U oobv);

method_t get_method(char const* buffer)
{
  // Remove '*' from the start
  if (*buffer == '*')
    ++buffer;

  // Remove 'bi' from the start
  if (buffer[0] == 'b' && buffer[1] == 'i')
    buffer += 2;

  switch (*buffer) {
    case 'n':
      return nearest;
    case 'l':
      return linear;
    case 'c':
      return cubic;
  }

  mexErrMsgTxt("Unsupported interpolation method");
  return nearest; // suppress warning.
}

template <class T>
bool try_cast(mlx_inputs& in, mlx_outputs& out);
  
template <class T, class Out_t>
bool try_cast_2(mlx_inputs& in, mlx_outputs& out);
  
void mlx_function(mlx_inputs& in, mlx_outputs& out)
{
  mlx_assert(in.nrhs >= 3 && in.nrhs <= 5);
  mlx_assert(out.nlhs <= 1);
  
//  mexPrintf("omp procs = %d\n", omp_get_num_procs());

  if (!(try_cast<mlx_logical>(in, out) ||
          try_cast<mlx_uint8>(in, out) ||
          try_cast<mlx_uint16>(in, out) ||
          try_cast<mlx_int16>(in, out) ||
          try_cast<mlx_single>(in, out) ||
          try_cast<mlx_double>(in, out)))
    mexErrMsgTxt("Unhandled input type (first argument)");
}

template <class T>
bool try_cast(mlx_inputs& in, mlx_outputs& out)
{
  // Check input type, fail silently if no match
  if (!mlx_isa<T>(in[0]))
    return false;

  // Check output type
  if (in.nrhs < 5)
    // It's double 
    mlx_assert((try_cast_2<T, mlx_double>(in, out)));
  else {
    if (!(try_cast_2<T, mlx_logical>(in, out) ||
            try_cast_2<T, mlx_uint8>(in, out) ||
            try_cast_2<T, mlx_uint16>(in, out) ||
            try_cast_2<T, mlx_single>(in, out) ||
            try_cast_2<T, mlx_double>(in, out)))
    mexErrMsgTxt("Onhandled output type (from oobv)");
  }
  return true;
} 
          
template <class T, class Out_t>
bool try_cast_2(mlx_inputs& in, mlx_outputs& out)
{
  // Get input, type should match already
  mlx_array<T> A(in[0]);
  mlx_assert(A.size.n <= 20);
  
  // Get the out of bounds value (oobv), should match Out_t as checked above.
  Out_t oobv;
  if (in.nrhs >= 5) {
    mlx_scalar<Out_t> oobv_array(in[4], mlx_array_nothrow);
    if (!oobv_array) return false;
    oobv = oobv_array[0];
  } else
    oobv = (Out_t)mxGetNaN();
  
  bool x = (bool)1.0;
  
  mlx_array<mlx_double> X(in[1]);
  mlx_array<mlx_double> Y(in[2]);
  
  mlx_assert(X.size == Y.size);
  mlx_assert(X.size.n <= 2);

  // Get the interpolation method
  method_t method = linear;
  if (in.nrhs > 3) {
    mlx_string s(in[3]);
    method = get_method(s.c_str());
  }
  
  int out_dim_buf[20];
  out_dim_buf[0] = X.size[0];
  out_dim_buf[1] = X.size[1];
  int nchannels = 1;
  for (int i = 2; i < A.size.n; i++) {
    out_dim_buf[i] = A.size[i];
    nchannels *= A.size[i];
  }
  mlx_dims out_dims(A.size.n, out_dim_buf);
  
  mlx_make_array<Out_t> ans(out_dims);

  switch (method) {
    case nearest:
      interp_nearest(ans, A, X, Y, oobv);
      break;
    case linear:
      interp_linear(ans, A, X, Y, oobv);
      break;
    case cubic:
      interp_cubic(ans, A, X, Y, oobv);
      break;
  }
  
  out[0] = ans;
  return true;
}

// Function for correct rounding
// Add these to use numeric_limits class
#include <limits>
using namespace std;
template<class U, class T> 
inline U saturate_cast(T val)
{
  if (numeric_limits<U>::is_integer && !numeric_limits<T>::is_integer) {
    if (numeric_limits<U>::is_signed)
      return val > 0 ? (val > (T)numeric_limits<U>::max() ? numeric_limits<U>::max() : static_cast<U>(val + 0.5)) : (val < (T)numeric_limits<U>::min() ? numeric_limits<U>::min() : static_cast<U>(val - 0.5));
    else
      return val > 0 ? (val > (T)numeric_limits<U>::max() ? numeric_limits<U>::max() : static_cast<U>(val + 0.5)) : 0;
  }
  return static_cast<U>(val);
}

// Nearest neighbour
template<class T, class Out_t>
void interp_nearest(mlx_array<Out_t>& B, mlx_array<T>& A, mlx_array<mlx_double>& X, mlx_array<mlx_double>& Y, Out_t oobv)
{
  int num_points = X.numel();

  int h = A.size[0];
  int w = A.size[1];
  int col = A.numel() / (w*h);

  int end = num_points * col;
  int step = h * w;

  double dw = (double)w + 0.5;
  double dh = (double)h + 0.5;
  
  // For each of the interpolation points
  int i, j, k;
#pragma omp parallel for if (num_points > 1000) num_threads(omp_get_num_procs()) default(shared) private(i,j,k)
  for (i = 0; i < num_points; i++) {
    double Xi = X[i];
    double Yi = Y[i];
    
    if (Xi >= 0.5 && Xi < dw && Yi >= 0.5 && Yi < dh) {
      // Find nearest neighbour
      k = h * int(Xi-0.5) + int(Yi-0.5);
      //mexPrintf("%d: %g %g -> %g\n", i, Xi, Yi, A[k]);
      for (j = i; j < end; j += num_points, k += step)
        B[j] = saturate_cast<Out_t, T>(A[k]);
    } else {
      // Out of bounds
      for (j = i; j < end; j += num_points)
        B[j] = oobv;
    }
  }
}

// Linear interpolation
template<class T, class Out_t>
void interp_linear(mlx_array<Out_t>& B, mlx_array<T>& A, mlx_array<mlx_double>& X, mlx_array<mlx_double>& Y, Out_t oobv)
{
  int num_points = X.numel();

  int h = A.size[0];
  int w = A.size[1];
  int col = A.numel() / (w*h);
  
  int end = num_points * col;
  int step = h * w;

  double dw = (double)w;
  double dh = (double)h;

  // For each of the interpolation points
  int i, j, k, x, y;
  double u, v, out;
#pragma omp parallel for if (num_points > 300) num_threads(omp_get_num_procs()) default(shared) private(i,j,k,u,v,x,y,out)
  for (i = 0; i < num_points; i++) {
    
    if (X[i] >= 1 && Y[i] >= 1) {
      if (X[i] < dw) {
        if (Y[i] < dh) {
          // Linearly interpolate
          x = (int)X[i];
          y = (int)Y[i];
          u = X[i] - x;
          v = Y[i] - y;
          k = h * (x - 1) + y - 1;
          for (j = i; j < end; j += num_points, k += step) {
            out = A[k] + (A[k+h] - A[k]) * u;
            out += ((A[k+1] - out) + (A[k+h+1] - A[k+1]) * u) * v;
            B[j] = saturate_cast<Out_t, double>(out);
          }
        } else if (Y[i] == dh) {
          // The Y coordinate is on the boundary
          // Avoid reading outside the buffer to avoid crashes
          // Linearly interpolate along X
          x = (int)X[i];
          u = X[i] - x;
          k = h * x - 1;
          for (j = i; j < end; j += num_points, k += step)
            B[j] = saturate_cast<Out_t, double>(A[k] + (A[k+h] - A[k]) * u);
        } else {
          // Out of bounds
          for (j = i; j < end; j += num_points)
            B[j] = oobv;
        }
      } else if (X[i] == dw) {
        if (Y[i] < dh) {
          // The X coordinate is on the boundary
          // Avoid reading outside the buffer to avoid crashes
          // Linearly interpolate along Y
          y = (int)Y[i];
          v = Y[i] - y;
          k = h * (w - 1) + y - 1;
          for (j = i; j < end; j += num_points, k += step)
            B[j] = saturate_cast<Out_t, double>(A[k] + (A[k+1] - A[k]) * v);
        } else if (Y[i] == dh) {
          // The X and Y coordinates are on the boundary
          // Avoid reading outside the buffer to avoid crashes
          // Output the last value in the array
          k = h * w - 1;
          for (j = i; j < end; j += num_points, k += step)
            B[j] = saturate_cast<Out_t, double>(A[k]);
        } else {
          // Out of bounds
          for (j = i; j < end; j += num_points)
            B[j] = oobv;
        }
      } else {
        // Out of bounds
        for (j = i; j < end; j += num_points)
          B[j] = oobv;
      }
    } else {
      // Out of bounds
      for (j = i; j < end; j += num_points)
        B[j] = oobv;
    }
  }
  return;
}

// Hermite cubic spline interpolation
template<class T, class Out_t>
void interp_cubic(mlx_array<Out_t>& B, mlx_array<T>& A, mlx_array<mlx_double>& X, mlx_array<mlx_double>& Y, Out_t oobv)
{
  int num_points = X.numel();

  int h = A.size[0];
  int w = A.size[1];
  int col = A.numel() / (w*h);
  mexPrintf("col = %d\n", col);
  int end = num_points * col;
  int step = h * w;

  double dw = (double)w - 1;
  double dh = (double)h - 1;

  // For each of the interpolation points
  int i, j, k, m, n, x, y;
  double a, b[4], c[4], u[3], v[3];
#pragma omp parallel for if (num_points > 100) num_threads(omp_get_num_procs()) default(shared) private(a,b,c,i,j,k,m,n,u,v,x,y)
  for (i = 0; i < num_points; i++) {
    
    if (X[i] >= 2 && X[i] < dw && Y[i] >= 2 && Y[i] < dh) {
      // Bicubicly interpolate
      x = (int)X[i];
      y = (int)Y[i];
      u[0] = X[i] - x;
      v[0] = Y[i] - y;
      u[1] = u[0] * u[0];
      v[1] = v[0] * v[0];
      u[2] = u[1] * u[0];
      v[2] = v[1] * v[0];
      k = h * (x - 2) + y - 2;
      for(j = i; j < end; j += num_points, k += step) {
        for (m = 0, n = k; m < 4; m++, n += h) {
          c[0] = (double)A[n+0];
          c[1] = (double)A[n+1];
          c[2] = (double)A[n+2];
          c[3] = (double)A[n+3];
          a = (c[3] + c[1]) - (c[2] + c[0]);
          b[m] = v[2] * a + v[1] * ((c[0] - c[1]) - a) + v[0] * (c[2] - c[0]) + c[1];
        }
        a = (b[3] + b[1]) - (b[2] + b[0]);
        B[j] = saturate_cast<Out_t, double>(u[2] * a + u[1] * ((b[0] - b[1]) - a) + u[0] * (b[2] - b[0]) + b[1]);
      }
    } else {
      // Out of bounds
      for (j = i; j < end; j += num_points)
        B[j] = oobv;
    }
  }
}
