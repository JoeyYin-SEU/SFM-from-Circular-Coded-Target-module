//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: minOrMax.h
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

#ifndef MINORMAX_H
#define MINORMAX_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
namespace internal {
float maximum(const float x[3]);

double maximum(const ::coder::array<double, 1U> &x);

void maximum(const ::coder::array<double, 2U> &x,
             ::coder::array<double, 2U> &ex);

void minimum(const ::coder::array<float, 1U> &x, float *ex, int *idx);

float minimum(const ::coder::array<float, 1U> &x);

} // namespace internal
} // namespace coder

#endif
//
// File trailer for minOrMax.h
//
// [EOF]
//
