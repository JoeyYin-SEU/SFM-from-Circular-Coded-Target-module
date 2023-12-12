//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: polyval.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

// Include Files
#include "polyval.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"

// Function Definitions
//
// Arguments    : coder::array<double, 2U> &in1
//                const coder::array<double, 2U> &in2
//                const double in3_data[]
//                int in4
// Return Type  : void
//
void binary_expand_op(coder::array<double, 2U> &in1,
                      const coder::array<double, 2U> &in2,
                      const double in3_data[], int in4)
{
  coder::array<double, 2U> b_in2;
  double in3;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  in3 = in3_data[in4 + 1];
  if (in1.size(1) == 1) {
    stride_0_1 = in2.size(1);
  } else {
    stride_0_1 = in1.size(1);
  }
  b_in2.set_size(1, stride_0_1);
  stride_0_1 = (in2.size(1) != 1);
  stride_1_1 = (in1.size(1) != 1);
  if (in1.size(1) == 1) {
    loop_ub = in2.size(1);
  } else {
    loop_ub = in1.size(1);
  }
  if (static_cast<int>(loop_ub < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      b_in2[i] = in2[i * stride_0_1] * in1[i * stride_1_1] + in3;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      b_in2[i] = in2[i * stride_0_1] * in1[i * stride_1_1] + in3;
    }
  }
  in1.set_size(1, b_in2.size(1));
  loop_ub = b_in2.size(1);
  if (static_cast<int>(b_in2.size(1) < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      in1[i] = b_in2[i];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      in1[i] = b_in2[i];
    }
  }
}

//
// File trailer for polyval.cpp
//
// [EOF]
//
