//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: sumprod.h
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

#ifndef SUMPROD_H
#define SUMPROD_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
void binary_expand_op(coder::array<int, 2U> &in1,
                      const coder::array<double, 2U> &in2,
                      const coder::array<double, 2U> &in3,
                      const coder::array<double, 2U> &in4);

void binary_expand_op(coder::array<int, 2U> &in1,
                      const coder::array<double, 2U> &in2,
                      const coder::array<double, 1U> &in3,
                      const coder::array<double, 2U> &in4);

#endif
//
// File trailer for sumprod.h
//
// [EOF]
//
