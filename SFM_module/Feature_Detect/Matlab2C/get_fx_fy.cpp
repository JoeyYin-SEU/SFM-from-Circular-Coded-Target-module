//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: get_fx_fy.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "get_fx_fy.h"
#include "combineVectorElements.h"
#include "get_chessborad_pixel_data.h"
#include "get_chessborad_pixel_initialize.h"
#include "imfilter.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>

// Function Definitions
//
// Arguments    : const coder::array<unsigned char, 2U> &b_I
//                double sigma
//                coder::array<float, 2U> &dx
//                coder::array<float, 2U> &dy
// Return Type  : void
//
void get_fx_fy(const coder::array<unsigned char, 2U> &b_I, double sigma,
               coder::array<float, 2U> &dx, coder::array<float, 2U> &dy)
{
  coder::array<double, 2U> b_derivGaussKernel;
  coder::array<double, 2U> derivGaussKernel;
  coder::array<double, 2U> x;
  coder::array<double, 1U> b_x;
  coder::array<int, 2U> r;
  coder::array<int, 2U> r1;
  coder::array<boolean_T, 2U> negVals;
  double a;
  double c;
  double filterExtent;
  double varargin_1;
  int loop_ub;
  int nx;
  if (!isInitialized_get_chessborad_pixel) {
    get_chessborad_pixel_initialize();
  }
  dy.set_size(b_I.size(0), b_I.size(1));
  loop_ub = b_I.size(0) * b_I.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int k{0}; k < loop_ub; k++) {
      dy[k] = static_cast<float>(b_I[k]) / 255.0F;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k < loop_ub; k++) {
      dy[k] = static_cast<float>(b_I[k]) / 255.0F;
    }
  }
  //  Create an even-length 1-D separable Derivative of Gaussian filter
  //  Determine filter length
  filterExtent = std::ceil(4.0 * sigma);
  if (std::isnan(-filterExtent) || std::isnan(filterExtent)) {
    x.set_size(1, 1);
    x[0] = rtNaN;
  } else if (filterExtent < -filterExtent) {
    x.set_size(x.size(0), 0);
  } else if ((std::isinf(-filterExtent) || std::isinf(filterExtent)) &&
             (-filterExtent == filterExtent)) {
    x.set_size(1, 1);
    x[0] = rtNaN;
  } else {
    c = -filterExtent;
    loop_ub = static_cast<int>(filterExtent - (-filterExtent));
    x.set_size(1, loop_ub + 1);
    if (static_cast<int>(loop_ub + 1 < 3200)) {
      for (int k{0}; k <= loop_ub; k++) {
        x[k] = -filterExtent + static_cast<double>(k);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int k = 0; k <= loop_ub; k++) {
        x[k] = c + static_cast<double>(k);
      }
    }
  }
  //  Create 1-D Gaussian Kernel
  a = 1.0 / (2.5066282746310002 * sigma);
  c = sigma * sigma;
  x.set_size(1, x.size(1));
  nx = x.size(1) - 1;
  loop_ub = x.size(1) - 1;
  if (static_cast<int>(x.size(1) < 3200)) {
    for (int k{0}; k <= nx; k++) {
      filterExtent = x[k];
      x[k] = filterExtent * filterExtent;
    }
  } else {
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(varargin_1)

    for (int k = 0; k <= loop_ub; k++) {
      varargin_1 = x[k];
      x[k] = varargin_1 * varargin_1;
    }
  }
  x.set_size(1, x.size(1));
  c *= 2.0;
  nx = x.size(1) - 1;
  loop_ub = x.size(1) - 1;
  if (static_cast<int>(x.size(1) < 3200)) {
    for (int k{0}; k <= nx; k++) {
      x[k] = -x[k] / c;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k <= loop_ub; k++) {
      x[k] = -x[k] / c;
    }
  }
  nx = x.size(1);
  if (static_cast<int>(x.size(1) < 3200)) {
    for (int k{0}; k < nx; k++) {
      x[k] = std::exp(x[k]);
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k < nx; k++) {
      x[k] = std::exp(x[k]);
    }
  }
  x.set_size(1, x.size(1));
  nx = x.size(1) - 1;
  loop_ub = x.size(1) - 1;
  if (static_cast<int>(x.size(1) < 3200)) {
    for (int k{0}; k <= nx; k++) {
      x[k] = a * x[k];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k <= loop_ub; k++) {
      x[k] = a * x[k];
    }
  }
  //  Normalize to ensure kernel sums to one
  filterExtent = coder::combineVectorElements(x);
  x.set_size(1, x.size(1));
  nx = x.size(1) - 1;
  loop_ub = x.size(1) - 1;
  if (static_cast<int>(x.size(1) < 3200)) {
    for (int k{0}; k <= nx; k++) {
      x[k] = x[k] / filterExtent;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k <= loop_ub; k++) {
      x[k] = x[k] / filterExtent;
    }
  }
  //  Create 1-D Derivative of Gaussian Kernel
  if (x.size(1) == 1) {
    derivGaussKernel.set_size(1, 1);
    derivGaussKernel[0] = 0.0;
  } else if (x.size(1) < 2) {
    derivGaussKernel.set_size(
        1, static_cast<int>(static_cast<signed char>(x.size(1))));
    loop_ub = static_cast<signed char>(x.size(1));
    for (nx = 0; nx < loop_ub; nx++) {
      derivGaussKernel[0] = 0.0;
    }
  } else {
    derivGaussKernel.set_size(1, x.size(1));
    derivGaussKernel[0] = x[1] - x[0];
    nx = x.size(1) - 1;
    if (static_cast<int>(x.size(1) - 2 < 3200)) {
      for (int k{2}; k <= nx; k++) {
        derivGaussKernel[k - 1] = (x[k] - x[k - 2]) / 2.0;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int k = 2; k <= nx; k++) {
        derivGaussKernel[k - 1] = (x[k] - x[k - 2]) / 2.0;
      }
    }
    derivGaussKernel[x.size(1) - 1] = x[x.size(1) - 1] - x[x.size(1) - 2];
  }
  //  Normalize to ensure kernel sums to zero
  negVals.set_size(1, derivGaussKernel.size(1));
  loop_ub = derivGaussKernel.size(1);
  if (static_cast<int>(derivGaussKernel.size(1) < 3200)) {
    for (int k{0}; k < loop_ub; k++) {
      negVals[k] = (derivGaussKernel[k] < 0.0);
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k < loop_ub; k++) {
      negVals[k] = (derivGaussKernel[k] < 0.0);
    }
  }
  loop_ub = derivGaussKernel.size(1) - 1;
  nx = 0;
  for (int i{0}; i <= loop_ub; i++) {
    if (derivGaussKernel[i] > 0.0) {
      nx++;
    }
  }
  r.set_size(1, nx);
  nx = 0;
  for (int i{0}; i <= loop_ub; i++) {
    if (derivGaussKernel[i] > 0.0) {
      r[nx] = i + 1;
      nx++;
    }
  }
  b_derivGaussKernel.set_size(1, r.size(1));
  loop_ub = r.size(1);
  if (static_cast<int>(r.size(1) < 3200)) {
    for (int k{0}; k < loop_ub; k++) {
      b_derivGaussKernel[k] = derivGaussKernel[r[k] - 1];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k < loop_ub; k++) {
      b_derivGaussKernel[k] = derivGaussKernel[r[k] - 1];
    }
  }
  filterExtent = coder::combineVectorElements(b_derivGaussKernel);
  b_derivGaussKernel.set_size(1, r.size(1));
  loop_ub = r.size(1);
  if (static_cast<int>(r.size(1) < 3200)) {
    for (int k{0}; k < loop_ub; k++) {
      b_derivGaussKernel[k] = derivGaussKernel[r[k] - 1] / filterExtent;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k < loop_ub; k++) {
      b_derivGaussKernel[k] = derivGaussKernel[r[k] - 1] / filterExtent;
    }
  }
  loop_ub = derivGaussKernel.size(1);
  nx = 0;
  for (int i{0}; i < loop_ub; i++) {
    if (derivGaussKernel[i] > 0.0) {
      derivGaussKernel[i] = b_derivGaussKernel[nx];
      nx++;
    }
  }
  loop_ub = negVals.size(1) - 1;
  nx = 0;
  for (int i{0}; i <= loop_ub; i++) {
    if (negVals[i]) {
      nx++;
    }
  }
  r1.set_size(1, nx);
  nx = 0;
  for (int i{0}; i <= loop_ub; i++) {
    if (negVals[i]) {
      r1[nx] = i + 1;
      nx++;
    }
  }
  b_derivGaussKernel.set_size(1, r1.size(1));
  loop_ub = r1.size(1);
  if (static_cast<int>(r1.size(1) < 3200)) {
    for (int k{0}; k < loop_ub; k++) {
      b_derivGaussKernel[k] = derivGaussKernel[r1[k] - 1];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k < loop_ub; k++) {
      b_derivGaussKernel[k] = derivGaussKernel[r1[k] - 1];
    }
  }
  filterExtent = std::abs(coder::combineVectorElements(b_derivGaussKernel));
  b_derivGaussKernel.set_size(1, r1.size(1));
  loop_ub = r1.size(1);
  if (static_cast<int>(r1.size(1) < 3200)) {
    for (int k{0}; k < loop_ub; k++) {
      b_derivGaussKernel[k] = derivGaussKernel[r1[k] - 1] / filterExtent;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k < loop_ub; k++) {
      b_derivGaussKernel[k] = derivGaussKernel[r1[k] - 1] / filterExtent;
    }
  }
  loop_ub = negVals.size(1);
  nx = 0;
  for (int i{0}; i < loop_ub; i++) {
    if (negVals[i]) {
      derivGaussKernel[i] = b_derivGaussKernel[nx];
      nx++;
    }
  }
  //  Compute smoothed numerical gradient of image I along x (horizontal)
  //  direction. GX corresponds to dG/dx, where G is the Gaussian Smoothed
  //  version of image I.
  dx.set_size(dy.size(0), dy.size(1));
  loop_ub = dy.size(0) * dy.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int k{0}; k < loop_ub; k++) {
      dx[k] = dy[k];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k < loop_ub; k++) {
      dx[k] = dy[k];
    }
  }
  b_x.set_size(x.size(1));
  loop_ub = x.size(1);
  if (static_cast<int>(x.size(1) < 3200)) {
    for (int k{0}; k < loop_ub; k++) {
      b_x[k] = x[k];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k < loop_ub; k++) {
      b_x[k] = x[k];
    }
  }
  coder::imfilter(dx, b_x);
  coder::imfilter(dx, derivGaussKernel);
  //  Compute smoothed numerical gradient of image I along y (vertical)
  //  direction. GY corresponds to dG/dy, where G is the Gaussian Smoothed
  //  version of image I.
  coder::imfilter(dy, x);
  b_x.set_size(derivGaussKernel.size(1));
  loop_ub = derivGaussKernel.size(1);
  if (static_cast<int>(derivGaussKernel.size(1) < 3200)) {
    for (int k{0}; k < loop_ub; k++) {
      b_x[k] = derivGaussKernel[k];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k < loop_ub; k++) {
      b_x[k] = derivGaussKernel[k];
    }
  }
  coder::imfilter(dy, b_x);
}

//
// File trailer for get_fx_fy.cpp
//
// [EOF]
//
