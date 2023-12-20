//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: squeeze.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "squeeze.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"

// Function Definitions
//
// Arguments    : const ::coder::array<double, 3U> &a
//                ::coder::array<double, 2U> &b
// Return Type  : void
//
namespace coder {
void b_squeeze(const ::coder::array<double, 3U> &a,
               ::coder::array<double, 2U> &b)
{
  int szb[2];
  int j;
  szb[0] = a.size(0);
  szb[1] = 1;
  if (a.size(2) != 1) {
    j = 0;
    if (a.size(0) != 1) {
      j = 1;
      szb[0] = a.size(0);
    }
    if (a.size(2) != 1) {
      szb[j] = a.size(2);
    }
  }
  b.set_size(szb[0], szb[1]);
  j = szb[0] * szb[1];
  if (static_cast<int>(j < 3200)) {
    for (int i{0}; i < j; i++) {
      b[i] = a[i];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < j; i++) {
      b[i] = a[i];
    }
  }
}

//
// Arguments    : const ::coder::array<double, 3U> &a
//                ::coder::array<double, 2U> &b
// Return Type  : void
//
void squeeze(const ::coder::array<double, 3U> &a, ::coder::array<double, 2U> &b)
{
  int szb[2];
  int j;
  szb[0] = 1;
  szb[1] = a.size(1);
  if (a.size(2) != 1) {
    j = 0;
    if (a.size(1) != 1) {
      j = 1;
      szb[0] = a.size(1);
    }
    if (a.size(2) != 1) {
      szb[j] = a.size(2);
    }
  }
  b.set_size(szb[0], szb[1]);
  j = szb[0] * szb[1];
  if (static_cast<int>(j < 3200)) {
    for (int i{0}; i < j; i++) {
      b[i] = a[i];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < j; i++) {
      b[i] = a[i];
    }
  }
}

} // namespace coder

//
// File trailer for squeeze.cpp
//
// [EOF]
//
