//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: polyfit.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "polyfit.h"
#include "qrsolve.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"
#include <algorithm>

// Function Definitions
//
// Arguments    : const ::coder::array<double, 1U> &x
//                const ::coder::array<double, 1U> &y
//                double n
//                double p_data[]
//                int p_size[2]
// Return Type  : void
//
namespace coder {
void polyfit(const ::coder::array<double, 1U> &x,
             const ::coder::array<double, 1U> &y, double n, double p_data[],
             int p_size[2])
{
  array<double, 2U> V;
  double p1_data[5];
  int p1_size;
  int rr;
  V.set_size(x.size(0), static_cast<int>(n + 1.0));
  if (x.size(0) != 0) {
    rr = x.size(0);
    for (int k{0}; k < rr; k++) {
      V[k + V.size(0) * (static_cast<int>(n + 1.0) - 1)] = 1.0;
    }
    rr = x.size(0);
    if (static_cast<int>(x.size(0) < 3200)) {
      for (int b_k{0}; b_k < rr; b_k++) {
        V[b_k + V.size(0) * (static_cast<int>(n) - 1)] = x[b_k];
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int b_k = 0; b_k < rr; b_k++) {
        V[b_k + V.size(0) * (static_cast<int>(n) - 1)] = x[b_k];
      }
    }
    rr = static_cast<int>(-((-1.0 - (n - 1.0)) + 1.0));
    p1_size = x.size(0);
    for (int j{0}; j < rr; j++) {
      double b_j;
      b_j = (n - 1.0) - static_cast<double>(j);
      for (int k{0}; k < p1_size; k++) {
        V[k + V.size(0) * (static_cast<int>(b_j) - 1)] =
            x[k] * V[k + V.size(0) * (static_cast<int>(b_j + 1.0) - 1)];
      }
    }
  }
  internal::qrsolve(V, y, p1_data, &p1_size, &rr);
  p_size[0] = 1;
  p_size[1] = p1_size;
  if (p1_size - 1 >= 0) {
    std::copy(&p1_data[0], &p1_data[p1_size], &p_data[0]);
  }
}

//
// Arguments    : const ::coder::array<double, 2U> &x
//                const ::coder::array<double, 2U> &y
//                double n
//                double p_data[]
//                int p_size[2]
// Return Type  : void
//
void polyfit(const ::coder::array<double, 2U> &x,
             const ::coder::array<double, 2U> &y, double n, double p_data[],
             int p_size[2])
{
  array<double, 2U> V;
  array<double, 1U> c_y;
  double p1_data[5];
  int b_y;
  int rr;
  V.set_size(x.size(1), static_cast<int>(n + 1.0));
  if (x.size(1) != 0) {
    rr = x.size(1);
    for (int k{0}; k < rr; k++) {
      V[k + V.size(0) * (static_cast<int>(n + 1.0) - 1)] = 1.0;
    }
    rr = x.size(1);
    if (static_cast<int>(x.size(1) < 3200)) {
      for (int b_k{0}; b_k < rr; b_k++) {
        V[b_k + V.size(0) * (static_cast<int>(n) - 1)] = x[b_k];
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int b_k = 0; b_k < rr; b_k++) {
        V[b_k + V.size(0) * (static_cast<int>(n) - 1)] = x[b_k];
      }
    }
    rr = static_cast<int>(-((-1.0 - (n - 1.0)) + 1.0));
    b_y = x.size(1);
    for (int j{0}; j < rr; j++) {
      double b_j;
      b_j = (n - 1.0) - static_cast<double>(j);
      for (int k{0}; k < b_y; k++) {
        V[k + V.size(0) * (static_cast<int>(b_j) - 1)] =
            x[k] * V[k + V.size(0) * (static_cast<int>(b_j + 1.0) - 1)];
      }
    }
  }
  b_y = y.size(1);
  c_y = y.reshape(b_y);
  internal::qrsolve(V, c_y, p1_data, &b_y, &rr);
  p_size[0] = 1;
  p_size[1] = b_y;
  if (b_y - 1 >= 0) {
    std::copy(&p1_data[0], &p1_data[b_y], &p_data[0]);
  }
}

} // namespace coder

//
// File trailer for polyfit.cpp
//
// [EOF]
//
