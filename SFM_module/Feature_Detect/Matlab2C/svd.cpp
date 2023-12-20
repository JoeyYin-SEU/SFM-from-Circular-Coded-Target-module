//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: svd.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "svd.h"
#include "rt_nonfinite.h"
#include "svd1.h"
#include "xzsvdc.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>
#include <cstring>

// Function Definitions
//
// Arguments    : const double A[841]
//                double U[841]
//                double S[841]
//                double V[841]
// Return Type  : void
//
namespace coder {
void svd(const double A[841], double U[841], double S[841], double V[841])
{
  double s[29];
  boolean_T p;
  p = true;
  for (int i{0}; i < 841; i++) {
    if (p) {
      double d;
      d = A[i];
      if (std::isinf(d) || std::isnan(d)) {
        p = false;
      }
    } else {
      p = false;
    }
  }
  if (p) {
    internal::b_svd(A, U, s, V);
  } else {
    for (int i{0}; i < 841; i++) {
      U[i] = rtNaN;
    }
    for (int i{0}; i < 29; i++) {
      s[i] = rtNaN;
    }
    for (int i{0}; i < 841; i++) {
      V[i] = rtNaN;
    }
  }
  std::memset(&S[0], 0, 841U * sizeof(double));
  for (int i{0}; i < 29; i++) {
    S[i + 29 * i] = s[i];
  }
}

//
// Arguments    : const ::coder::array<double, 2U> &A
//                double *U
//                ::coder::array<double, 2U> &S
//                ::coder::array<double, 2U> &V
// Return Type  : void
//
void svd(const ::coder::array<double, 2U> &A, double *U,
         ::coder::array<double, 2U> &S, ::coder::array<double, 2U> &V)
{
  array<double, 2U> V1;
  array<double, 2U> r;
  double s_data;
  double x;
  int m;
  int nx;
  boolean_T p;
  nx = A.size(1);
  p = true;
  for (int k{0}; k < nx; k++) {
    if (p) {
      x = A[k];
      if (std::isinf(x) || std::isnan(x)) {
        p = false;
      }
    } else {
      p = false;
    }
  }
  if (p) {
    if (A.size(1) == 0) {
      V.set_size(0, 0);
      *U = 1.0;
      m = 0;
    } else {
      internal::reflapack::xzsvdc(A, U, (double *)&s_data, &m, V);
    }
  } else {
    unsigned int unnamed_idx_1;
    unnamed_idx_1 = static_cast<unsigned int>(A.size(1));
    if (A.size(1) == 0) {
      m = A.size(1);
      V1.set_size(A.size(1), A.size(1));
      nx = A.size(1) * A.size(1);
      if (static_cast<int>(nx < 3200)) {
        for (int i{0}; i < nx; i++) {
          V1[i] = 0.0;
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int i = 0; i < nx; i++) {
          V1[i] = 0.0;
        }
      }
      if (static_cast<int>(unnamed_idx_1) > 0) {
        for (int k{0}; k < m; k++) {
          V1[k + V1.size(0) * k] = 1.0;
        }
      }
      m = 0;
    } else {
      r.set_size(1, A.size(1));
      nx = A.size(1);
      if (static_cast<int>(A.size(1) < 3200)) {
        for (int i{0}; i < nx; i++) {
          r[i] = 0.0;
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int i = 0; i < nx; i++) {
          r[i] = 0.0;
        }
      }
      internal::reflapack::xzsvdc(r, &x, (double *)&s_data, &m, V1);
    }
    *U = rtNaN;
    for (nx = 0; nx < m; nx++) {
      s_data = rtNaN;
    }
    V.set_size(V1.size(0), V1.size(1));
    nx = V1.size(0) * V1.size(1);
    if (static_cast<int>(nx < 3200)) {
      for (int i{0}; i < nx; i++) {
        V[i] = rtNaN;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < nx; i++) {
        V[i] = rtNaN;
      }
    }
  }
  S.set_size(1, V.size(1));
  nx = V.size(1);
  if (static_cast<int>(V.size(1) < 3200)) {
    for (int i{0}; i < nx; i++) {
      S[i] = 0.0;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < nx; i++) {
      S[i] = 0.0;
    }
  }
  if (m - 1 >= 0) {
    S[0] = s_data;
  }
}

//
// Arguments    : const ::coder::array<double, 1U> &A
//                ::coder::array<double, 2U> &U
//                ::coder::array<double, 1U> &S
//                double *V
// Return Type  : void
//
void svd(const ::coder::array<double, 1U> &A, ::coder::array<double, 2U> &U,
         ::coder::array<double, 1U> &S, double *V)
{
  array<double, 2U> U1;
  array<double, 1U> r;
  double V1;
  double s_data;
  int nx;
  boolean_T p;
  nx = A.size(0);
  p = true;
  for (int k{0}; k < nx; k++) {
    if ((!p) || (std::isinf(A[k]) || std::isnan(A[k]))) {
      p = false;
    }
  }
  if (p) {
    internal::b_svd(A, U, (double *)&s_data, &nx, V);
  } else {
    r.set_size(A.size(0));
    nx = A.size(0);
    if (static_cast<int>(A.size(0) < 3200)) {
      for (int i{0}; i < nx; i++) {
        r[i] = 0.0;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < nx; i++) {
        r[i] = 0.0;
      }
    }
    internal::b_svd(r, U1, (double *)&s_data, &nx, &V1);
    U.set_size(U1.size(0), U1.size(1));
    nx = U1.size(0) * U1.size(1);
    if (static_cast<int>(nx < 3200)) {
      for (int i{0}; i < nx; i++) {
        U[i] = rtNaN;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < nx; i++) {
        U[i] = rtNaN;
      }
    }
    s_data = rtNaN;
    *V = rtNaN;
  }
  S.set_size(U.size(1));
  nx = U.size(1);
  if (static_cast<int>(U.size(1) < 3200)) {
    for (int i{0}; i < nx; i++) {
      S[i] = 0.0;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < nx; i++) {
      S[i] = 0.0;
    }
  }
  S[0] = s_data;
}

} // namespace coder

//
// File trailer for svd.cpp
//
// [EOF]
//
