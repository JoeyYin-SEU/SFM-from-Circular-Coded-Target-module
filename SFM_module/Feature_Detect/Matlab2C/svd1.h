//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: svd1.h
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

#ifndef SVD1_H
#define SVD1_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
namespace internal {
void b_svd(const double A[841], double U[841], double s[29], double V[841]);

void b_svd(const ::coder::array<double, 1U> &A, ::coder::array<double, 2U> &U,
           double s_data[], int *s_size, double *V);

} // namespace internal
} // namespace coder

#endif
//
// File trailer for svd1.h
//
// [EOF]
//
