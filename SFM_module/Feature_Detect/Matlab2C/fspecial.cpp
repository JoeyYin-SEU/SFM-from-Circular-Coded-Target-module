//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: fspecial.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

// Include Files
#include "fspecial.h"
#include "rt_nonfinite.h"
#include "omp.h"

// Function Definitions
//
// Arguments    : double in1_data[]
//                int in1_size[2]
//                const double in2_data[]
//                const int in2_size[2]
//                double in3
// Return Type  : void
//
void binary_expand_op(double in1_data[], int in1_size[2],
                      const double in2_data[], const int in2_size[2],
                      double in3)
{
  double b_in1_data[225];
  double d;
  int b_in1_size[2];
  int aux_0_1;
  int aux_1_1;
  int c_loop_ub;
  int i3;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  d = 2.0 * in3 * in3;
  if (in2_size[0] == 1) {
    b_in1_size[0] = in1_size[0];
  } else {
    b_in1_size[0] = in2_size[0];
  }
  if (in2_size[1] == 1) {
    b_in1_size[1] = in1_size[1];
  } else {
    b_in1_size[1] = in2_size[1];
  }
  stride_0_0 = (in1_size[0] != 1);
  stride_0_1 = (in1_size[1] != 1);
  stride_1_0 = (in2_size[0] != 1);
  stride_1_1 = (in2_size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (in2_size[1] == 1) {
    loop_ub = in1_size[1];
  } else {
    loop_ub = in2_size[1];
  }
  for (int i{0}; i < loop_ub; i++) {
    int b_loop_ub;
    int i1;
    i1 = in2_size[0];
    if (i1 == 1) {
      b_loop_ub = in1_size[0];
    } else {
      b_loop_ub = i1;
    }
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      double in1_tmp;
      int in1_data_tmp;
      in1_tmp = in2_data[i1 * stride_1_0 + in2_size[0] * aux_1_1];
      in1_data_tmp = in1_size[0] * aux_0_1;
      b_in1_data[i1 + b_in1_size[0] * i] =
          -(in1_data[i1 * stride_0_0 + in1_data_tmp] *
                in1_data[i1 * stride_0_0 + in1_data_tmp] +
            in1_tmp * in1_tmp) /
          d;
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  in1_size[0] = b_in1_size[0];
  in1_size[1] = b_in1_size[1];
  loop_ub = b_in1_size[1];
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(i3, c_loop_ub)

  for (int i2 = 0; i2 < loop_ub; i2++) {
    c_loop_ub = b_in1_size[0];
    for (i3 = 0; i3 < c_loop_ub; i3++) {
      in1_data[i3 + in1_size[0] * i2] = b_in1_data[i3 + b_in1_size[0] * i2];
    }
  }
}

//
// File trailer for fspecial.cpp
//
// [EOF]
//
