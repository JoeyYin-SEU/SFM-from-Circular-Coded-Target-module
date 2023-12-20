//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: qrsolve.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "qrsolve.h"
#include "rt_nonfinite.h"
#include "xgeqp3.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>
#include <cstring>

// Function Definitions
//
// Arguments    : const ::coder::array<double, 2U> &A
//                const ::coder::array<double, 1U> &B
//                double Y_data[]
//                int *Y_size
//                int *rankA
// Return Type  : void
//
namespace coder {
namespace internal {
void qrsolve(const ::coder::array<double, 2U> &A,
             const ::coder::array<double, 1U> &B, double Y_data[], int *Y_size,
             int *rankA)
{
  array<double, 2U> b_A;
  array<double, 1U> b_B;
  double tau_data[5];
  double tol;
  int jpvt_data[5];
  int jpvt_size[2];
  int assumedRank;
  int m;
  int maxmn;
  int minmn;
  int u1;
  b_A.set_size(A.size(0), A.size(1));
  minmn = A.size(0) * A.size(1);
  if (static_cast<int>(minmn < 3200)) {
    for (int i{0}; i < minmn; i++) {
      b_A[i] = A[i];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < minmn; i++) {
      b_A[i] = A[i];
    }
  }
  lapack::xgeqp3(b_A, tau_data, &minmn, jpvt_data, jpvt_size);
  *rankA = 0;
  if (b_A.size(0) < b_A.size(1)) {
    minmn = b_A.size(0);
    maxmn = b_A.size(1);
  } else {
    minmn = b_A.size(1);
    maxmn = b_A.size(0);
  }
  if (minmn > 0) {
    tol = std::fmin(1.4901161193847656E-8,
                    2.2204460492503131E-15 * static_cast<double>(maxmn)) *
          std::abs(b_A[0]);
    while ((*rankA < minmn) &&
           (!(std::abs(b_A[*rankA + b_A.size(0) * *rankA]) <= tol))) {
      (*rankA)++;
    }
  }
  assumedRank = 0;
  maxmn = b_A.size(0);
  minmn = b_A.size(1);
  if (maxmn <= minmn) {
    minmn = maxmn;
  }
  if (minmn > 0) {
    for (maxmn = 0; maxmn < minmn; maxmn++) {
      if (b_A[maxmn + b_A.size(0) * maxmn] != 0.0) {
        assumedRank++;
      }
    }
  }
  b_B.set_size(B.size(0));
  minmn = B.size(0);
  if (static_cast<int>(B.size(0) < 3200)) {
    for (int i{0}; i < minmn; i++) {
      b_B[i] = B[i];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < minmn; i++) {
      b_B[i] = B[i];
    }
  }
  *Y_size = b_A.size(1);
  minmn = b_A.size(1);
  if (minmn - 1 >= 0) {
    std::memset(&Y_data[0], 0,
                static_cast<unsigned int>(minmn) * sizeof(double));
  }
  m = b_A.size(0);
  maxmn = b_A.size(0);
  u1 = b_A.size(1);
  if (maxmn <= u1) {
    u1 = maxmn;
  }
  for (maxmn = 0; maxmn < u1; maxmn++) {
    if (tau_data[maxmn] != 0.0) {
      tol = b_B[maxmn];
      minmn = maxmn + 2;
      for (int b_i{minmn}; b_i <= m; b_i++) {
        tol += b_A[(b_i + b_A.size(0) * maxmn) - 1] * b_B[b_i - 1];
      }
      tol *= tau_data[maxmn];
      if (tol != 0.0) {
        b_B[maxmn] = b_B[maxmn] - tol;
        for (int b_i{minmn}; b_i <= m; b_i++) {
          b_B[b_i - 1] =
              b_B[b_i - 1] - b_A[(b_i + b_A.size(0) * maxmn) - 1] * tol;
        }
      }
    }
  }
  for (int b_i{0}; b_i < assumedRank; b_i++) {
    Y_data[jpvt_data[b_i] - 1] = b_B[b_i];
  }
  for (maxmn = assumedRank; maxmn >= 1; maxmn--) {
    minmn = jpvt_data[maxmn - 1];
    Y_data[minmn - 1] /= b_A[(maxmn + b_A.size(0) * (maxmn - 1)) - 1];
    for (int b_i{0}; b_i <= maxmn - 2; b_i++) {
      Y_data[jpvt_data[b_i] - 1] -= Y_data[jpvt_data[maxmn - 1] - 1] *
                                    b_A[b_i + b_A.size(0) * (maxmn - 1)];
    }
  }
}

} // namespace internal
} // namespace coder

//
// File trailer for qrsolve.cpp
//
// [EOF]
//
