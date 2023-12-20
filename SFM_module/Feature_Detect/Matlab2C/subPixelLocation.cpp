//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: subPixelLocation.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "subPixelLocation.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <algorithm>
#include <cmath>

// Variable Definitions
static double X[150];

// Function Declarations
static void binary_expand_op(float in1_data[], int in1_size[2],
                             const coder::array<double, 2U> &in2, int in3,
                             float in4, float in5);

// Function Definitions
//
// Arguments    : float in1_data[]
//                int in1_size[2]
//                const coder::array<double, 2U> &in2
//                int in3
//                float in4
//                float in5
// Return Type  : void
//
static void binary_expand_op(float in1_data[], int in1_size[2],
                             const coder::array<double, 2U> &in2, int in3,
                             float in4, float in5)
{
  in1_size[0] = 1;
  in1_size[1] = 2;
  in1_data[0] = static_cast<float>(in2[in3]) + in4;
  in1_data[1] =
      static_cast<float>(in2[in3 + in2.size(0) * (in2.size(1) != 1)]) + in5;
}

//
// Arguments    : const ::coder::array<float, 2U> &metric
//                ::coder::array<double, 2U> &loc
// Return Type  : void
//
namespace coder {
namespace vision {
namespace internal {
namespace calibration {
namespace checkerboard {
void subPixelLocation(const ::coder::array<float, 2U> &metric,
                      ::coder::array<double, 2U> &loc)
{
  array<float, 2U> b_metric;
  array<boolean_T, 1U> x;
  int subPixelLoc_size[2];
  int i;
  i = loc.size(0);
  for (int id{0}; id < i; id++) {
    float subPixelLoc_data[2];
    int b_loc;
    int i1;
    int loop_ub;
    boolean_T exitg1;
    boolean_T y;
    b_loc = loc.size(1);
    x.set_size(b_loc);
    for (i1 = 0; i1 < b_loc; i1++) {
      x[i1] = (loc[id + loc.size(0) * i1] < 3.0);
    }
    y = false;
    b_loc = 1;
    exitg1 = false;
    while ((!exitg1) && (b_loc <= x.size(0))) {
      if (x[b_loc - 1]) {
        y = true;
        exitg1 = true;
      } else {
        b_loc++;
      }
    }
    if (y || (loc[id] > (static_cast<double>(metric.size(1)) - 2.0) - 1.0) ||
        (loc[id + loc.size(0)] >
         (static_cast<double>(metric.size(0)) - 2.0) - 1.0)) {
      subPixelLoc_size[0] = 1;
      loop_ub = loc.size(1);
      subPixelLoc_size[1] = loop_ub;
      for (i1 = 0; i1 < loop_ub; i1++) {
        subPixelLoc_data[i1] = static_cast<float>(loc[id + loc.size(0) * i1]);
      }
    } else {
      float beta[6];
      float b_x;
      float b_y;
      int b_loop_ub;
      int i2;
      int i3;
      if (loc[id + loc.size(0)] - 2.0 > loc[id + loc.size(0)] + 2.0) {
        i1 = 0;
        i2 = 0;
      } else {
        i1 = static_cast<int>(loc[id + loc.size(0)] - 2.0) - 1;
        i2 = static_cast<int>(loc[id + loc.size(0)] + 2.0);
      }
      if (loc[id] - 2.0 > loc[id] + 2.0) {
        b_loc = 0;
        i3 = 0;
      } else {
        b_loc = static_cast<int>(loc[id] - 2.0) - 1;
        i3 = static_cast<int>(loc[id] + 2.0);
      }
      loop_ub = i2 - i1;
      b_loop_ub = i3 - b_loc;
      b_metric.set_size(loop_ub, b_loop_ub);
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        for (i3 = 0; i3 < loop_ub; i3++) {
          b_metric[i3 + b_metric.size(0) * i2] =
              metric[(i1 + i3) + metric.size(0) * (b_loc + i2)];
        }
      }
      for (i1 = 0; i1 < 6; i1++) {
        b_y = 0.0F;
        for (i2 = 0; i2 < 25; i2++) {
          b_y += static_cast<float>(X[i1 + 6 * i2]) * b_metric[i2];
        }
        beta[i1] = b_y;
      }
      b_y = 4.0F * beta[0] * beta[1] - beta[4] * beta[4];
      b_x = -(2.0F * beta[1] * beta[2] - beta[3] * beta[4]) / b_y;
      b_y = -(2.0F * beta[0] * beta[3] - beta[2] * beta[4]) / b_y;
      if (std::isinf(b_x) || std::isnan(b_x) || (std::abs(b_x) > 2.0F) ||
          (std::isinf(b_y) || std::isnan(b_y)) || (std::abs(b_y) > 2.0F)) {
        b_x = 0.0F;
        b_y = 0.0F;
      }
      if (loc.size(1) == 2) {
        subPixelLoc_size[0] = 1;
        subPixelLoc_size[1] = 2;
        subPixelLoc_data[0] = static_cast<float>(loc[id]) + b_x;
        subPixelLoc_data[1] = static_cast<float>(loc[id + loc.size(0)]) + b_y;
      } else {
        binary_expand_op(subPixelLoc_data, subPixelLoc_size, loc, id, b_x, b_y);
      }
    }
    loop_ub = subPixelLoc_size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      loc[id + loc.size(0) * i1] = subPixelLoc_data[i1];
    }
  }
}

//
// Arguments    : void
// Return Type  : void
//
} // namespace checkerboard
} // namespace calibration
} // namespace internal
} // namespace vision
} // namespace coder
void subPixelLocationImpl_init()
{
  static const double dv[150]{0.028571428571428574,
                              0.028571428571428574,
                              -0.04,
                              -0.04,
                              0.04,
                              -0.0742857142857143,
                              0.028571428571428574,
                              -0.014285714285714285,
                              -0.04,
                              -0.02,
                              0.02,
                              0.01142857142857142,
                              0.028571428571428574,
                              -0.028571428571428571,
                              -0.04,
                              0.0,
                              0.0,
                              0.039999999999999994,
                              0.028571428571428574,
                              -0.014285714285714285,
                              -0.04,
                              0.02,
                              -0.02,
                              0.01142857142857142,
                              0.028571428571428574,
                              0.028571428571428574,
                              -0.04,
                              0.04,
                              -0.04,
                              -0.0742857142857143,
                              -0.014285714285714287,
                              0.028571428571428571,
                              -0.02,
                              -0.04,
                              0.02,
                              0.011428571428571429,
                              -0.014285714285714285,
                              -0.014285714285714284,
                              -0.02,
                              -0.02,
                              0.01,
                              0.097142857142857142,
                              -0.01428571428571429,
                              -0.028571428571428574,
                              -0.02,
                              0.0,
                              0.0,
                              0.12571428571428572,
                              -0.014285714285714285,
                              -0.014285714285714284,
                              -0.02,
                              0.02,
                              -0.01,
                              0.097142857142857142,
                              -0.014285714285714287,
                              0.028571428571428571,
                              -0.02,
                              0.04,
                              -0.02,
                              0.011428571428571429,
                              -0.028571428571428574,
                              0.028571428571428571,
                              0.0,
                              -0.04,
                              0.0,
                              0.040000000000000008,
                              -0.028571428571428574,
                              -0.014285714285714287,
                              0.0,
                              -0.02,
                              0.0,
                              0.12571428571428572,
                              -0.028571428571428574,
                              -0.028571428571428574,
                              0.0,
                              0.0,
                              0.0,
                              0.1542857142857143,
                              -0.028571428571428574,
                              -0.014285714285714287,
                              0.0,
                              0.02,
                              0.0,
                              0.12571428571428572,
                              -0.028571428571428574,
                              0.028571428571428571,
                              0.0,
                              0.04,
                              0.0,
                              0.040000000000000008,
                              -0.014285714285714287,
                              0.028571428571428571,
                              0.02,
                              -0.04,
                              -0.02,
                              0.011428571428571429,
                              -0.014285714285714285,
                              -0.014285714285714284,
                              0.02,
                              -0.02,
                              -0.01,
                              0.097142857142857142,
                              -0.01428571428571429,
                              -0.028571428571428574,
                              0.02,
                              0.0,
                              0.0,
                              0.12571428571428572,
                              -0.014285714285714285,
                              -0.014285714285714284,
                              0.02,
                              0.02,
                              0.01,
                              0.097142857142857142,
                              -0.014285714285714287,
                              0.028571428571428571,
                              0.02,
                              0.04,
                              0.02,
                              0.011428571428571429,
                              0.028571428571428574,
                              0.028571428571428574,
                              0.04,
                              -0.04,
                              -0.04,
                              -0.0742857142857143,
                              0.028571428571428574,
                              -0.014285714285714285,
                              0.04,
                              -0.02,
                              -0.02,
                              0.01142857142857142,
                              0.028571428571428574,
                              -0.028571428571428571,
                              0.04,
                              0.0,
                              0.0,
                              0.039999999999999994,
                              0.028571428571428574,
                              -0.014285714285714285,
                              0.04,
                              0.02,
                              0.02,
                              0.01142857142857142,
                              0.028571428571428574,
                              0.028571428571428574,
                              0.04,
                              0.04,
                              0.04,
                              -0.0742857142857143};
  std::copy(&dv[0], &dv[150], &X[0]);
}

//
// File trailer for subPixelLocation.cpp
//
// [EOF]
//
