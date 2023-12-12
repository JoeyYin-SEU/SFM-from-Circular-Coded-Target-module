//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: svd.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

// Include Files
#include "svd.h"
#include "rt_nonfinite.h"
#include "svd1.h"
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

} // namespace coder

//
// File trailer for svd.cpp
//
// [EOF]
//
