//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: svd.h
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

#ifndef SVD_H
#define SVD_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
void svd(const double A[841], double U[841], double S[841], double V[841]);

void svd(const ::coder::array<double, 2U> &A, double *U,
         ::coder::array<double, 2U> &S, ::coder::array<double, 2U> &V);

void svd(const ::coder::array<double, 1U> &A, ::coder::array<double, 2U> &U,
         ::coder::array<double, 1U> &S, double *V);

} // namespace coder

#endif
//
// File trailer for svd.h
//
// [EOF]
//
