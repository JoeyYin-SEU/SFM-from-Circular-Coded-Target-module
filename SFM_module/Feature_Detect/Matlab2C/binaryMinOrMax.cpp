//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: binaryMinOrMax.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "binaryMinOrMax.h"
#include "minOrMax.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>

// Function Definitions
//
// Arguments    : float in1
//                const coder::array<double, 1U> &in2
//                const coder::array<double, 1U> &in3
// Return Type  : float
//
float binary_expand_op(float in1, const coder::array<double, 1U> &in2,
                       const coder::array<double, 1U> &in3)
{
  coder::array<double, 1U> b_in2;
  float out1;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  if (in3.size(0) == 1) {
    stride_0_0 = in2.size(0);
  } else {
    stride_0_0 = in3.size(0);
  }
  b_in2.set_size(stride_0_0);
  stride_0_0 = (in2.size(0) != 1);
  stride_1_0 = (in3.size(0) != 1);
  if (in3.size(0) == 1) {
    loop_ub = in2.size(0);
  } else {
    loop_ub = in3.size(0);
  }
  if (static_cast<int>(loop_ub < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      b_in2[i] = in2[i * stride_0_0] / in3[i * stride_1_0];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      b_in2[i] = in2[i * stride_0_0] / in3[i * stride_1_0];
    }
  }
  out1 = std::fmax(in1, static_cast<float>(coder::internal::maximum(b_in2)));
  return out1;
}

//
// File trailer for binaryMinOrMax.cpp
//
// [EOF]
//
