#include <mex.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Assert macros
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

// ----------------------------------------------------------------------------

// This is a (thankfully not too) smart pointer to an mxArray.
template <class T>
struct mlx_cast {
  mxClassID matlab_class_id;
  bool ok;
  T* data;
  mwSize rows;
  mwSize cols;
  mwSize number_of_dims;
  mwSize const *dims;
  
  mlx_cast(mxArray const* a) {
    matlab_class_id = mlx_class_id((T*)0);
    ok = (a) && (mxGetClassID(a) == matlab_class_id);
    if (!ok) return;
    
    // This is the correct type, set up the data
    data = (T*)mxGetData(a);
    rows = mxGetM(a);
    cols = mxGetN(a);
    number_of_dims = mxGetNumberOfDimensions(a);
    dims =  mxGetDimensions(a);
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
};
