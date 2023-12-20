//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: unsafeSxfun.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "unsafeSxfun.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"

// Function Definitions
//
// Arguments    : coder::array<float, 2U> &in1
//                const coder::array<double, 2U> &in3
//                const coder::array<float, 2U> &in4
// Return Type  : void
//
void binary_expand_op(coder::array<float, 2U> &in1,
                      const coder::array<double, 2U> &in3,
                      const coder::array<float, 2U> &in4)
{
  coder::array<float, 2U> b_in3;
  float b_varargin_1;
  int aux_1_1;
  int i;
  int i3;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_1_1;
  if (in4.size(0) == 1) {
    i = in3.size(0);
  } else {
    i = in4.size(0);
  }
  b_in3.set_size(i, 2);
  stride_0_0 = (in3.size(0) != 1);
  stride_1_0 = (in4.size(0) != 1);
  stride_1_1 = (in4.size(1) != 1);
  aux_1_1 = 0;
  if (in4.size(0) == 1) {
    loop_ub = in3.size(0);
  } else {
    loop_ub = in4.size(0);
  }
  for (i = 0; i < 2; i++) {
    for (int i1{0}; i1 < loop_ub; i1++) {
      b_in3[i1 + b_in3.size(0) * i] =
          static_cast<float>(in3[i1 * stride_0_0 + in3.size(0) * i]) -
          in4[i1 * stride_1_0 + in4.size(0) * aux_1_1];
    }
    aux_1_1 += stride_1_1;
  }
  in1.set_size(b_in3.size(0), 2);
  loop_ub = b_in3.size(0);
  if (static_cast<int>((b_in3.size(0) << 1) < 3200)) {
    for (int i2{0}; i2 < 2; i2++) {
      for (i3 = 0; i3 < loop_ub; i3++) {
        float varargin_1;
        varargin_1 = b_in3[i3 + b_in3.size(0) * i2];
        in1[i3 + in1.size(0) * i2] = varargin_1 * varargin_1;
      }
    }
  } else {
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(i3, b_varargin_1)

    for (int i2 = 0; i2 < 2; i2++) {
      for (i3 = 0; i3 < loop_ub; i3++) {
        b_varargin_1 = b_in3[i3 + b_in3.size(0) * i2];
        in1[i3 + in1.size(0) * i2] = b_varargin_1 * b_varargin_1;
      }
    }
  }
}

//
// File trailer for unsafeSxfun.cpp
//
// [EOF]
//
