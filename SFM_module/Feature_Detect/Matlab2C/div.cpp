//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: div.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "div.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"

// Function Definitions
//
// Arguments    : coder::array<float, 1U> &in1
//                const coder::array<float, 1U> &in2
//                float in3
// Return Type  : void
//
void binary_expand_op(coder::array<float, 1U> &in1,
                      const coder::array<float, 1U> &in2, float in3)
{
  coder::array<float, 1U> b_in1;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  if (in2.size(0) == 1) {
    stride_0_0 = in1.size(0);
  } else {
    stride_0_0 = in2.size(0);
  }
  b_in1.set_size(stride_0_0);
  stride_0_0 = (in1.size(0) != 1);
  stride_1_0 = (in2.size(0) != 1);
  if (in2.size(0) == 1) {
    loop_ub = in1.size(0);
  } else {
    loop_ub = in2.size(0);
  }
  if (static_cast<int>(loop_ub < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      b_in1[i] = in1[i * stride_0_0] / (in2[i * stride_1_0] * in3);
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      b_in1[i] = in1[i * stride_0_0] / (in2[i * stride_1_0] * in3);
    }
  }
  in1.set_size(b_in1.size(0));
  loop_ub = b_in1.size(0);
  if (static_cast<int>(b_in1.size(0) < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      in1[i] = b_in1[i];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      in1[i] = b_in1[i];
    }
  }
}

//
// File trailer for div.cpp
//
// [EOF]
//
