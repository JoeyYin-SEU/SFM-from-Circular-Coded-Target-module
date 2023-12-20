//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: sumprod.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "sumprod.h"
#include "combineVectorElements.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

// Function Definitions
//
// Arguments    : coder::array<int, 2U> &in1
//                const coder::array<double, 2U> &in2
//                const coder::array<double, 2U> &in3
//                const coder::array<double, 2U> &in4
// Return Type  : void
//
void binary_expand_op(coder::array<int, 2U> &in1,
                      const coder::array<double, 2U> &in2,
                      const coder::array<double, 2U> &in3,
                      const coder::array<double, 2U> &in4)
{
  coder::array<boolean_T, 2U> b_in2;
  int aux_0_1;
  int aux_1_1;
  int b_in3;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  b_in3 = in3.size(1);
  if (in4.size(1) == 1) {
    i = in2.size(1);
  } else {
    i = in4.size(1);
  }
  b_in2.set_size(in2.size(0), i);
  stride_0_1 = (in2.size(1) != 1);
  stride_1_1 = (in4.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (in4.size(1) == 1) {
    loop_ub = in2.size(1);
  } else {
    loop_ub = in4.size(1);
  }
  for (i = 0; i < loop_ub; i++) {
    int b_loop_ub;
    b_loop_ub = in2.size(0);
    for (int i1{0}; i1 < b_loop_ub; i1++) {
      b_in2[i1 + b_in2.size(0) * i] =
          (in2[i1 + in2.size(0) * aux_0_1] >
           static_cast<double>(b_in3) * in4[aux_1_1] * 2.2204460492503131E-16);
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  coder::b_combineVectorElements(b_in2, in1);
}

//
// Arguments    : coder::array<int, 2U> &in1
//                const coder::array<double, 2U> &in2
//                const coder::array<double, 1U> &in3
//                const coder::array<double, 2U> &in4
// Return Type  : void
//
void binary_expand_op(coder::array<int, 2U> &in1,
                      const coder::array<double, 2U> &in2,
                      const coder::array<double, 1U> &in3,
                      const coder::array<double, 2U> &in4)
{
  coder::array<boolean_T, 2U> b_in2;
  int aux_0_1;
  int aux_1_1;
  int b_in3;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  b_in3 = in3.size(0);
  if (in4.size(1) == 1) {
    i = in2.size(1);
  } else {
    i = in4.size(1);
  }
  b_in2.set_size(in2.size(0), i);
  stride_0_1 = (in2.size(1) != 1);
  stride_1_1 = (in4.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (in4.size(1) == 1) {
    loop_ub = in2.size(1);
  } else {
    loop_ub = in4.size(1);
  }
  for (i = 0; i < loop_ub; i++) {
    int b_loop_ub;
    b_loop_ub = in2.size(0);
    for (int i1{0}; i1 < b_loop_ub; i1++) {
      b_in2[i1 + b_in2.size(0) * i] =
          (in2[i1 + in2.size(0) * aux_0_1] >
           static_cast<double>(b_in3) * in4[aux_1_1] * 2.2204460492503131E-16);
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  coder::b_combineVectorElements(b_in2, in1);
}

//
// File trailer for sumprod.cpp
//
// [EOF]
//
