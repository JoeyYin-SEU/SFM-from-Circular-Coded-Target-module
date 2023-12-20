//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: rot90.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "rot90.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"

// Function Definitions
//
// Arguments    : const ::coder::array<double, 2U> &A
//                ::coder::array<double, 2U> &B
// Return Type  : void
//
namespace coder {
void b_rot90(const ::coder::array<double, 2U> &A, ::coder::array<double, 2U> &B)
{
  int i;
  int m;
  int n;
  m = A.size(0);
  n = A.size(1);
  B.set_size(A.size(0), A.size(1));
  if (A.size(0) * A.size(1) >= 8192) {
    int ub_loop;
    ub_loop = A.size(1) - 1;
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32) private(i)

    for (int j = 0; j <= ub_loop; j++) {
      for (i = 0; i < m; i++) {
        B[i + B.size(0) * j] = A[((m - i) + A.size(0) * ((n - j) - 1)) - 1];
      }
    }
  } else if (static_cast<int>(A.size(0) * A.size(1) < 3200)) {
    for (int j{0}; j < n; j++) {
      for (i = 0; i < m; i++) {
        B[i + B.size(0) * j] = A[((m - i) + A.size(0) * ((n - j) - 1)) - 1];
      }
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32) private(i)

    for (int j = 0; j < n; j++) {
      for (i = 0; i < m; i++) {
        B[i + B.size(0) * j] = A[((m - i) + A.size(0) * ((n - j) - 1)) - 1];
      }
    }
  }
}

//
// Arguments    : const ::coder::array<double, 2U> &A
//                ::coder::array<double, 2U> &B
// Return Type  : void
//
void rot90(const ::coder::array<double, 2U> &A, ::coder::array<double, 2U> &B)
{
  int j;
  int m;
  int n;
  m = A.size(0);
  n = A.size(1);
  B.set_size(A.size(1), A.size(0));
  if (A.size(0) * A.size(1) >= 8192) {
    int ub_loop;
    ub_loop = A.size(1) - 1;
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32) private(j)

    for (int i = 0; i <= ub_loop; i++) {
      for (j = 0; j < m; j++) {
        B[i + B.size(0) * j] = A[j + A.size(0) * ((n - i) - 1)];
      }
    }
  } else if (static_cast<int>(A.size(0) * A.size(1) < 3200)) {
    for (int i{0}; i < n; i++) {
      for (j = 0; j < m; j++) {
        B[i + B.size(0) * j] = A[j + A.size(0) * ((n - i) - 1)];
      }
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32) private(j)

    for (int i = 0; i < n; i++) {
      for (j = 0; j < m; j++) {
        B[i + B.size(0) * j] = A[j + A.size(0) * ((n - i) - 1)];
      }
    }
  }
}

} // namespace coder

//
// File trailer for rot90.cpp
//
// [EOF]
//
