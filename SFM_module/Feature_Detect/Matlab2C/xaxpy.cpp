//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xaxpy.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

// Include Files
#include "xaxpy.h"
#include "rt_nonfinite.h"
#include "omp.h"

// Function Definitions
//
// Arguments    : int n
//                double a
//                const double x[29]
//                int ix0
//                double y[841]
//                int iy0
// Return Type  : void
//
namespace coder {
namespace internal {
namespace blas {
void b_xaxpy(int n, double a, const double x[29], int ix0, double y[841],
             int iy0)
{
  if (!(a == 0.0)) {
    int i;
    i = n - 1;
    for (int k{0}; k <= i; k++) {
      int i1;
      i1 = (iy0 + k) - 1;
      y[i1] += a * x[(ix0 + k) - 1];
    }
  }
}

//
// Arguments    : int n
//                double a
//                const double x[841]
//                int ix0
//                double y[29]
//                int iy0
// Return Type  : void
//
void xaxpy(int n, double a, const double x[841], int ix0, double y[29], int iy0)
{
  int i1;
  if (!(a == 0.0)) {
    int i;
    int ix;
    int iy;
    ix = ix0 - 1;
    iy = iy0 - 1;
    i = n - 1;
    if (static_cast<int>(n < 3200)) {
      for (int k{0}; k <= i; k++) {
        ix = (iy0 + k) - 1;
        y[ix] += a * x[(ix0 + k) - 1];
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32) private(i1)

      for (int k = 0; k <= i; k++) {
        i1 = iy + k;
        y[i1] += a * x[ix + k];
      }
    }
  }
}

//
// Arguments    : int n
//                double a
//                int ix0
//                double y[841]
//                int iy0
// Return Type  : void
//
void xaxpy(int n, double a, int ix0, double y[841], int iy0)
{
  if (!(a == 0.0)) {
    int i;
    i = n - 1;
    for (int k{0}; k <= i; k++) {
      int i1;
      i1 = (iy0 + k) - 1;
      y[i1] += a * y[(ix0 + k) - 1];
    }
  }
}

} // namespace blas
} // namespace internal
} // namespace coder

//
// File trailer for xaxpy.cpp
//
// [EOF]
//
