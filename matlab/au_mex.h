/*
 *  AU_MEX.H -- Simplify MEX file creation
 *
 * Mex files are a great way to speed up operations, but can be a 
 * little arcane to code.   This header file contains a number of routines 
 * to simplify them.
 *
 * The most important class is mlx_cast<T> which behaves like a pointer to
 * an mxArray.  To write a fully error-checked mex file which adds two 
 * double arrays, write this:
 *
 * 
 * */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <mex.h>

// Assert macro.
// This is not disabled by optimization -- if you want to make an assert
// in a tight loop, wrap it in #ifdef MLX_DEBUG
#define mlx_assert(expr) if (expr) 0; else mexErrMsgTxt("Assertion failed: " ## #expr)


// Declare C types for matlab types, e.g. mlx_uint32 or mlx_single
// Declare mlx_class_id(T*), mapping from declared types, e.g. mlx_int32, or int32_T, to classid
#define DECLARE_MEX_CLASS(ID, matlab_type, ctype) \
  typedef ctype mlx_ ## matlab_type; \
  inline mxClassID mlx_class_id(mlx_ ## matlab_type *) { return ID; }
  
DECLARE_MEX_CLASS(mxINT8_CLASS, int8, int8_T)
DECLARE_MEX_CLASS(mxUINT8_CLASS, uint8, uint8_T)
DECLARE_MEX_CLASS(mxINT16_CLASS, int16, int16_T)
DECLARE_MEX_CLASS(mxUINT16_CLASS, uint16, uint16_T)
DECLARE_MEX_CLASS(mxINT32_CLASS, int32, int32_T)
DECLARE_MEX_CLASS(mxUINT32_CLASS, uint32, uint32_T)
DECLARE_MEX_CLASS(mxINT64_CLASS, int64, int64_T)
DECLARE_MEX_CLASS(mxUINT64_CLASS, uint64, uint64_T)
DECLARE_MEX_CLASS(mxSINGLE_CLASS, single, float)
DECLARE_MEX_CLASS(mxDOUBLE_CLASS, double, double)
DECLARE_MEX_CLASS(mxLOGICAL_CLASS, logical, bool) // ????

template <class T>
struct mlx_isa {
  bool it_is;
  mlx_isa(mxArray const* a) {
    it_is = (mxGetClassID(a) == mlx_class_id((T*)0));
  }
  operator bool () const { return it_is; }
};



// ----------------------------------------------------------------------------

template <class T>
struct mlx_numeric_traits {
};

template <>
struct mlx_numeric_traits<unsigned int> {
  typedef unsigned int wide_t;
};

template <>
struct mlx_numeric_traits<mlx_uint8> {
  typedef unsigned int wide_t;
};

template <>
struct mlx_numeric_traits<mlx_single> {
  typedef mlx_double wide_t;
};

template <>
struct mlx_numeric_traits<mlx_double> {
  typedef mlx_double wide_t;
};

// ----------------------------------------------------------------------------

template <class T>
struct mlx_traits {
};

template <>
struct mlx_traits<mlx_uint8> {
  typedef mlx_uint8 T;
  typedef mlx_numeric_traits<T>::wide_t wide_t;
  static T colour_max() { return 255; }
  static T narrowing_divide_by_max(wide_t val) { return T(val / 255); }
};

template <>
struct mlx_traits<mlx_single> {
  typedef mlx_single T;
  typedef mlx_numeric_traits<T>::wide_t wide_t;

  static T colour_max() { return 1.0f; }
  static T narrowing_divide_by_max(wide_t t) { return T(t); }
};

template <>
struct mlx_traits<mlx_double> {
  typedef mlx_double T;
  typedef mlx_numeric_traits<T>::wide_t wide_t;

  static T colour_max() { return 1.0; }
  static T narrowing_divide_by_max(wide_t t) { return T(t); }
};

template <class T>
bool assert_equal(T*,T*) {};

template <class T1, class T2>
struct mlp_pairwise_traits {
  typedef typename mlx_numeric_traits<T1>::wide_t wide1_t;
  typedef typename mlx_numeric_traits<wide1_t>::wide_t wide_t;
  
  typedef typename mlx_numeric_traits<
    typename mlx_numeric_traits<T2>::wide_t
  >::wide_t wide2_t;
  void f() { assert_equal((wide_t*)0, (wide2_t*)0); }
};

// ----------------------------------------------------------------------------

template <class T1, class T2>
inline
typename mlp_pairwise_traits<T1,T2>::wide_t
widening_product(T1 a, T2 b)
{
  typedef typename mlp_pairwise_traits<T1,T2>::wide_t out_t;
  return out_t(out_t(a) * out_t(b));
}

inline unsigned int saturating_add(unsigned int a, unsigned int b)
{
  return a + b;
}


inline double saturating_add(double a, double b)
{
  return a + b;
}

// A struct to hold the mx dimensions.
// This does not need to take ownership of the pointer,
// as it is guaranteed to stay alive as long as the matrix does.
struct mlx_dims {
  mwSize n;
  mwSize const *dims;
  
  mlx_dims():n(0), dims(0) {}
  mlx_dims(mwSize n_, mwSize const* dims_):n(n_), dims(dims_) {}
//  mlx_dims(int n_, int* dims_):n(n_), dims(dims_) {}
//  mlx_dims(size_t n_, size_t const* dims_):n(n_), dims(dims_) {}
  mlx_dims(const mlx_dims& that):n(that.n), dims(that.dims) {}
  
};

// Check dims for equality
bool operator==(mlx_dims const& a, mlx_dims const& b)
{
    if (a.n != b.n) return false;
    for(int i = 0; i < a.n; ++i)
        if (a.dims[i] != b.dims[i])
            return false;
    
    return true;
}

// ----------------------------------------------------------------------------

// Pass to mlx_cast to indicate that type mismatch is non-fatal.
static const bool mlx_cast_nothrow = false;

// This is a (thankfully not too) smart pointer to an mxArray.
template <class T>
struct mlx_cast {
  mxArray const* mx_array;
  mxClassID matlab_class_id;
  bool ok;
  T* data;
  mwSize rows;
  mwSize cols;
  mlx_dims size;
  
  // Take mxArray pointer 'a', and check that its contents
  // correspond to the template type parameter T.
  // If so, take local copies of its size and data pointer.
  mlx_cast(mxArray const* a, bool throw_on_type_mismatch = true) 
  {
    mx_array = a;
    matlab_class_id = mlx_class_id((T*)0);
    ok = (a) && (mxGetClassID(a) == matlab_class_id);
    if (!ok) 
        if (throw_on_type_mismatch && a)
            mexErrMsgIdAndTxt("awful:bad_cast", "Unexpected datatype [%s]", mxGetClassName(a));
        else
            return; // Return silently, caller can check flag;
    
    // This is the correct type, set up the data
    data = (T*)mxGetData(a);
    rows = (mwSize)mxGetM(a);
    cols = (mwSize)mxGetN(a);
    size = mlx_dims(mxGetNumberOfDimensions(a), mxGetDimensions(a));
    mlx_assert(mxGetPi(a) == 0); // Don't handle complex yet.
  }
  
  // 1-based put
  void put1(mwIndex r, mwIndex c, const T& t) {
    data[r-1 + (c-1)*rows] = t;
  }
  // 1-based get
  T const& get1(mwIndex r, mwIndex c) const {
    return data[r-1 + (c-1)*rows];
  }

  // 0-based put
  void put0(mwIndex r, mwIndex c, const T& t) {
    data[r + c*rows] = t;
  }
  // 0-based get
  T const& get0(mwIndex r, mwIndex c) const {
    return data[r + c*rows];
  }

  // 0-based []
  T const& operator[](mwIndex ind) const {
    return data[ind];
  }

  // 0-based []
  T & operator[](mwIndex ind) {
    return data[ind];
  }

  // Check if the type cast in the constructor succeeded.
  operator bool() const { return ok; }
  
  // Return number of elements
  mwSize numel() const {
      if (size.n == 0)
          return 0;
      
      mwSize n = size.dims[0];
      for(mwSize i = 1; i < size.n; ++i)
          n *= size.dims[i];

      return n;
  }

};

// Create an mxArray of given size, with 
// the matlab class id corresponding to C++ type T
template <class T>
struct mlx_make_array : mlx_cast<T>
{
    // Construct from (rows, cols)
    mlx_make_array(mwSize rows, mwSize cols):
        base_t(0)
    {
        mwSize dims[2] = {rows, cols};
        mlx_dims sz(2, dims);
        create(sz);
    }

    // Construct from size array
    mlx_make_array(mlx_dims const& sz):
        base_t(0)
    {
        create(sz);
    }

private:
    typedef mlx_cast<T> base_t;
    
    void create(mlx_dims const& size_)
    {
        size = size_;
        int *odims = (int *)mxMalloc(sizeof(int)*size.n);
        for(int i=0; i<size.n; i++)
            odims[i] = size.dims[i];
        mx_array = mxCreateNumericArray(size.n, odims, 
                matlab_class_id, mxREAL);

        data = (T*)mxGetData(mx_array);
        rows = (mwSize)mxGetM(mx_array);
        cols = (mwSize)mxGetN(mx_array);
  }
};

// Class to collect input arguments, and apply error checking to their
// access.
struct mlx_inputs {
    int nrhs;
    mxArray const** prhs;
    
    // Construct from nrhs, prhs arguments to mexFunction
    mlx_inputs(int nrhs, mxArray const* prhs[]):nrhs(nrhs), prhs(prhs)
    {
    }
    
    mxArray const* operator[](int i) {
        mlx_assert(i >= 0);
        if (i >= nrhs)
            mexErrMsgIdAndTxt("awful:nin", "Expected at least %d input arguments", i+1);
        return prhs[i];
    }
};

// Helper class to allow
// out[i] = x;
// where x is an mlx_cast<> object;
struct mlx_output {
    mxArray** array_ptr;
    
    // Construct from pointer to mxArray*
    mlx_output(mxArray** array_ptr = 0):array_ptr(array_ptr) {}
    
    // Assign from mlx_cast<T>
    template <class T>
    mlx_output& operator=(mlx_cast<T>& that) {
        *array_ptr = const_cast<mxArray*>(that.mx_array);
        return *this;
    }
    
    // Assign from mxArray*
    template <class T>
    mlx_output& operator=(mxArray* that) {
        *array_ptr = that;
        return *this;
    }
};

// Class to collect output arguments, and apply error checking to their
// access.
struct mlx_outputs {
    int nlhs;
    mxArray** plhs;
    
    // Construct from nlhs, plhs arguments to mexFunction
    mlx_outputs(int nlhs, mxArray * plhs[]):nlhs(nlhs), plhs(plhs)
    {
    }
    // Access to ith output argument, with range checking
    // Special case i==0, as one may assign to that even if nlhs==0
    mlx_output operator[](int i)
    {
        mlx_assert(i >= 0);
        if (i > 0 && i >= nlhs)
            mexErrMsgIdAndTxt("awful:nout", "Expected at least %d output arguments", i+1);
        return mlx_output(&plhs[i]);
    }
};
