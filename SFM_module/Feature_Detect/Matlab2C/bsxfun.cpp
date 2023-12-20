//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: bsxfun.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "bsxfun.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"

// Function Definitions
//
// Arguments    : const ::coder::array<float, 3U> &a
//                const ::coder::array<float, 2U> &b
//                ::coder::array<float, 3U> &c
// Return Type  : void
//
namespace coder {
void bsxfun(const ::coder::array<float, 3U> &a,
            const ::coder::array<float, 2U> &b, ::coder::array<float, 3U> &c)
{
  int acoef;
  int b_acoef;
  int b_k;
  int bcoef;
  int c_k;
  int i;
  int i1;
  int varargin_3;
  int varargin_5;
  acoef = b.size(0);
  bcoef = a.size(0);
  if (acoef <= bcoef) {
    bcoef = acoef;
  }
  acoef = b.size(1);
  b_acoef = a.size(1);
  if (acoef <= b_acoef) {
    b_acoef = acoef;
  }
  if (b.size(0) == 1) {
    bcoef = a.size(0);
  } else if (a.size(0) == 1) {
    bcoef = b.size(0);
  } else if (a.size(0) == b.size(0)) {
    bcoef = a.size(0);
  }
  if (b.size(1) == 1) {
    b_acoef = a.size(1);
  } else if (a.size(1) == 1) {
    b_acoef = b.size(1);
  } else if (a.size(1) == b.size(1)) {
    b_acoef = a.size(1);
  }
  c.set_size(bcoef, b_acoef, 2);
  acoef = b.size(0);
  bcoef = a.size(0);
  if (acoef <= bcoef) {
    bcoef = acoef;
  }
  acoef = b.size(1);
  b_acoef = a.size(1);
  if (acoef <= b_acoef) {
    b_acoef = acoef;
  }
  if (b.size(0) == 1) {
    bcoef = a.size(0);
  } else if (a.size(0) == 1) {
    bcoef = b.size(0);
  } else if (a.size(0) == b.size(0)) {
    bcoef = a.size(0);
  }
  if (b.size(1) == 1) {
    b_acoef = a.size(1);
  } else if (a.size(1) == 1) {
    b_acoef = b.size(1);
  } else if (a.size(1) == b.size(1)) {
    b_acoef = a.size(1);
  }
  if ((bcoef != 0) && (b_acoef != 0)) {
    int b_bcoef;
    acoef = (a.size(1) != 1);
    bcoef = (b.size(1) != 1);
    b_acoef = (a.size(0) != 1);
    b_bcoef = (b.size(0) != 1);
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(b_k, i, varargin_5,     \
                                                       varargin_3, c_k, i1)

    for (int k = 0; k < 2; k++) {
      i1 = c.size(1) - 1;
      for (c_k = 0; c_k <= i1; c_k++) {
        varargin_3 = acoef * c_k;
        varargin_5 = bcoef * c_k;
        i = c.size(0) - 1;
        for (b_k = 0; b_k <= i; b_k++) {
          c[(b_k + c.size(0) * c_k) + c.size(0) * c.size(1) * k] =
              a[(b_acoef * b_k + a.size(0) * varargin_3) +
                a.size(0) * a.size(1) * k] -
              b[b_bcoef * b_k + b.size(0) * varargin_5];
        }
      }
    }
  }
}

//
// Arguments    : const ::coder::array<float, 2U> &a
//                const ::coder::array<double, 2U> &b
//                ::coder::array<float, 2U> &c
// Return Type  : void
//
void bsxfun(const ::coder::array<float, 2U> &a,
            const ::coder::array<double, 2U> &b, ::coder::array<float, 2U> &c)
{
  int acoef;
  int b_k;
  int csz_idx_1;
  int i;
  int u0;
  int varargin_2;
  int varargin_3;
  u0 = b.size(1);
  acoef = a.size(1);
  if (u0 <= acoef) {
    acoef = u0;
  }
  if (b.size(1) == 1) {
    csz_idx_1 = a.size(1);
  } else if (a.size(1) == 1) {
    csz_idx_1 = b.size(1);
  } else if (a.size(1) == b.size(1)) {
    csz_idx_1 = a.size(1);
  } else {
    csz_idx_1 = acoef;
  }
  u0 = b.size(1);
  acoef = a.size(1);
  if (u0 <= acoef) {
    acoef = u0;
  }
  if (b.size(1) == 1) {
    acoef = a.size(1);
  } else if (a.size(1) == 1) {
    acoef = b.size(1);
  } else if (a.size(1) == b.size(1)) {
    acoef = a.size(1);
  }
  c.set_size(a.size(0), acoef);
  if ((a.size(0) != 0) && (csz_idx_1 != 0)) {
    int b_acoef;
    int bcoef;
    b_acoef = (a.size(1) != 1);
    bcoef = (b.size(1) != 1);
    u0 = csz_idx_1 - 1;
    acoef = (a.size(0) != 1);
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads()                         \
                               : 32) private(b_k, i, varargin_3, varargin_2)

    for (int k = 0; k <= u0; k++) {
      varargin_2 = b_acoef * k;
      varargin_3 = bcoef * k;
      i = c.size(0) - 1;
      for (b_k = 0; b_k <= i; b_k++) {
        c[b_k + c.size(0) * k] = a[acoef * b_k + a.size(0) * varargin_2] -
                                 static_cast<float>(b[varargin_3]);
      }
    }
  }
}

} // namespace coder

//
// File trailer for bsxfun.cpp
//
// [EOF]
//
