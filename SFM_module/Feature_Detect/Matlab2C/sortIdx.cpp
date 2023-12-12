//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: sortIdx.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

// Include Files
#include "sortIdx.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>

// Function Declarations
namespace coder {
namespace internal {
static void merge(::coder::array<int, 1U> &idx, ::coder::array<float, 1U> &x,
                  int offset, int np, int nq, ::coder::array<int, 1U> &iwork,
                  ::coder::array<float, 1U> &xwork);

static void merge(::coder::array<int, 1U> &idx, ::coder::array<int, 1U> &x,
                  int offset, int np, int nq, ::coder::array<int, 1U> &iwork,
                  ::coder::array<int, 1U> &xwork);

static void merge(::coder::array<int, 1U> &idx, ::coder::array<double, 1U> &x,
                  int offset, int np, int nq, ::coder::array<int, 1U> &iwork,
                  ::coder::array<double, 1U> &xwork);

} // namespace internal
} // namespace coder

// Function Definitions
//
// Arguments    : ::coder::array<int, 1U> &idx
//                ::coder::array<float, 1U> &x
//                int offset
//                int np
//                int nq
//                ::coder::array<int, 1U> &iwork
//                ::coder::array<float, 1U> &xwork
// Return Type  : void
//
namespace coder {
namespace internal {
static void merge(::coder::array<int, 1U> &idx, ::coder::array<float, 1U> &x,
                  int offset, int np, int nq, ::coder::array<int, 1U> &iwork,
                  ::coder::array<float, 1U> &xwork)
{
  int i1;
  if (nq != 0) {
    int i;
    int iout;
    int n;
    int p;
    n = np + nq;
    if (static_cast<int>(n < 3200)) {
      for (int j{0}; j < n; j++) {
        i = offset + j;
        iwork[j] = idx[i];
        xwork[j] = x[i];
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32) private(i1)

      for (int j = 0; j < n; j++) {
        i1 = offset + j;
        iwork[j] = idx[i1];
        xwork[j] = x[i1];
      }
    }
    p = 0;
    n = np;
    i = np + nq;
    iout = offset - 1;
    int exitg1;
    do {
      exitg1 = 0;
      iout++;
      if (xwork[p] >= xwork[n]) {
        idx[iout] = iwork[p];
        x[iout] = xwork[p];
        if (p + 1 < np) {
          p++;
        } else {
          exitg1 = 1;
        }
      } else {
        idx[iout] = iwork[n];
        x[iout] = xwork[n];
        if (n + 1 < i) {
          n++;
        } else {
          n = iout - p;
          if (static_cast<int>(np - p < 3200)) {
            for (int j{p + 1}; j <= np; j++) {
              i = n + j;
              idx[i] = iwork[j - 1];
              x[i] = xwork[j - 1];
            }
          } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32) private(i1)

            for (int j = p + 1; j <= np; j++) {
              i1 = n + j;
              idx[i1] = iwork[j - 1];
              x[i1] = xwork[j - 1];
            }
          }
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

//
// Arguments    : ::coder::array<int, 1U> &idx
//                ::coder::array<int, 1U> &x
//                int offset
//                int np
//                int nq
//                ::coder::array<int, 1U> &iwork
//                ::coder::array<int, 1U> &xwork
// Return Type  : void
//
static void merge(::coder::array<int, 1U> &idx, ::coder::array<int, 1U> &x,
                  int offset, int np, int nq, ::coder::array<int, 1U> &iwork,
                  ::coder::array<int, 1U> &xwork)
{
  int i1;
  if (nq != 0) {
    int i;
    int iout;
    int n;
    int p;
    n = np + nq;
    if (static_cast<int>(n < 3200)) {
      for (int j{0}; j < n; j++) {
        i = offset + j;
        iwork[j] = idx[i];
        xwork[j] = x[i];
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32) private(i1)

      for (int j = 0; j < n; j++) {
        i1 = offset + j;
        iwork[j] = idx[i1];
        xwork[j] = x[i1];
      }
    }
    p = 0;
    n = np;
    i = np + nq;
    iout = offset - 1;
    int exitg1;
    do {
      exitg1 = 0;
      iout++;
      if (xwork[p] <= xwork[n]) {
        idx[iout] = iwork[p];
        x[iout] = xwork[p];
        if (p + 1 < np) {
          p++;
        } else {
          exitg1 = 1;
        }
      } else {
        idx[iout] = iwork[n];
        x[iout] = xwork[n];
        if (n + 1 < i) {
          n++;
        } else {
          n = iout - p;
          if (static_cast<int>(np - p < 3200)) {
            for (int j{p + 1}; j <= np; j++) {
              i = n + j;
              idx[i] = iwork[j - 1];
              x[i] = xwork[j - 1];
            }
          } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32) private(i1)

            for (int j = p + 1; j <= np; j++) {
              i1 = n + j;
              idx[i1] = iwork[j - 1];
              x[i1] = xwork[j - 1];
            }
          }
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

//
// Arguments    : ::coder::array<int, 1U> &idx
//                ::coder::array<double, 1U> &x
//                int offset
//                int np
//                int nq
//                ::coder::array<int, 1U> &iwork
//                ::coder::array<double, 1U> &xwork
// Return Type  : void
//
static void merge(::coder::array<int, 1U> &idx, ::coder::array<double, 1U> &x,
                  int offset, int np, int nq, ::coder::array<int, 1U> &iwork,
                  ::coder::array<double, 1U> &xwork)
{
  int i1;
  if (nq != 0) {
    int i;
    int iout;
    int n;
    int p;
    n = np + nq;
    if (static_cast<int>(n < 3200)) {
      for (int j{0}; j < n; j++) {
        i = offset + j;
        iwork[j] = idx[i];
        xwork[j] = x[i];
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32) private(i1)

      for (int j = 0; j < n; j++) {
        i1 = offset + j;
        iwork[j] = idx[i1];
        xwork[j] = x[i1];
      }
    }
    p = 0;
    n = np;
    i = np + nq;
    iout = offset - 1;
    int exitg1;
    do {
      exitg1 = 0;
      iout++;
      if (xwork[p] <= xwork[n]) {
        idx[iout] = iwork[p];
        x[iout] = xwork[p];
        if (p + 1 < np) {
          p++;
        } else {
          exitg1 = 1;
        }
      } else {
        idx[iout] = iwork[n];
        x[iout] = xwork[n];
        if (n + 1 < i) {
          n++;
        } else {
          n = iout - p;
          if (static_cast<int>(np - p < 3200)) {
            for (int j{p + 1}; j <= np; j++) {
              i = n + j;
              idx[i] = iwork[j - 1];
              x[i] = xwork[j - 1];
            }
          } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32) private(i1)

            for (int j = p + 1; j <= np; j++) {
              i1 = n + j;
              idx[i1] = iwork[j - 1];
              x[i1] = xwork[j - 1];
            }
          }
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

//
// Arguments    : ::coder::array<int, 1U> &idx
//                ::coder::array<float, 1U> &x
//                int offset
//                int n
//                int preSortLevel
//                ::coder::array<int, 1U> &iwork
//                ::coder::array<float, 1U> &xwork
// Return Type  : void
//
void merge_block(::coder::array<int, 1U> &idx, ::coder::array<float, 1U> &x,
                 int offset, int n, int preSortLevel,
                 ::coder::array<int, 1U> &iwork,
                 ::coder::array<float, 1U> &xwork)
{
  int bLen;
  int nPairs;
  nPairs = n >> preSortLevel;
  bLen = 1 << preSortLevel;
  while (nPairs > 1) {
    int nTail;
    int tailOffset;
    if ((nPairs & 1) != 0) {
      nPairs--;
      tailOffset = bLen * nPairs;
      nTail = n - tailOffset;
      if (nTail > bLen) {
        merge(idx, x, offset + tailOffset, bLen, nTail - bLen, iwork, xwork);
      }
    }
    tailOffset = bLen << 1;
    nPairs >>= 1;
    for (nTail = 0; nTail < nPairs; nTail++) {
      merge(idx, x, offset + nTail * tailOffset, bLen, bLen, iwork, xwork);
    }
    bLen = tailOffset;
  }
  if (n > bLen) {
    merge(idx, x, offset, bLen, n - bLen, iwork, xwork);
  }
}

//
// Arguments    : ::coder::array<int, 1U> &idx
//                ::coder::array<int, 1U> &x
//                int offset
//                int n
//                int preSortLevel
//                ::coder::array<int, 1U> &iwork
//                ::coder::array<int, 1U> &xwork
// Return Type  : void
//
void merge_block(::coder::array<int, 1U> &idx, ::coder::array<int, 1U> &x,
                 int offset, int n, int preSortLevel,
                 ::coder::array<int, 1U> &iwork, ::coder::array<int, 1U> &xwork)
{
  int bLen;
  int nPairs;
  nPairs = n >> preSortLevel;
  bLen = 1 << preSortLevel;
  while (nPairs > 1) {
    int nTail;
    int tailOffset;
    if ((nPairs & 1) != 0) {
      nPairs--;
      tailOffset = bLen * nPairs;
      nTail = n - tailOffset;
      if (nTail > bLen) {
        merge(idx, x, offset + tailOffset, bLen, nTail - bLen, iwork, xwork);
      }
    }
    tailOffset = bLen << 1;
    nPairs >>= 1;
    for (nTail = 0; nTail < nPairs; nTail++) {
      merge(idx, x, offset + nTail * tailOffset, bLen, bLen, iwork, xwork);
    }
    bLen = tailOffset;
  }
  if (n > bLen) {
    merge(idx, x, offset, bLen, n - bLen, iwork, xwork);
  }
}

//
// Arguments    : ::coder::array<int, 1U> &idx
//                ::coder::array<double, 1U> &x
//                int offset
//                int n
//                int preSortLevel
//                ::coder::array<int, 1U> &iwork
//                ::coder::array<double, 1U> &xwork
// Return Type  : void
//
void merge_block(::coder::array<int, 1U> &idx, ::coder::array<double, 1U> &x,
                 int offset, int n, int preSortLevel,
                 ::coder::array<int, 1U> &iwork,
                 ::coder::array<double, 1U> &xwork)
{
  int bLen;
  int nPairs;
  nPairs = n >> preSortLevel;
  bLen = 1 << preSortLevel;
  while (nPairs > 1) {
    int nTail;
    int tailOffset;
    if ((nPairs & 1) != 0) {
      nPairs--;
      tailOffset = bLen * nPairs;
      nTail = n - tailOffset;
      if (nTail > bLen) {
        merge(idx, x, offset + tailOffset, bLen, nTail - bLen, iwork, xwork);
      }
    }
    tailOffset = bLen << 1;
    nPairs >>= 1;
    for (nTail = 0; nTail < nPairs; nTail++) {
      merge(idx, x, offset + nTail * tailOffset, bLen, bLen, iwork, xwork);
    }
    bLen = tailOffset;
  }
  if (n > bLen) {
    merge(idx, x, offset, bLen, n - bLen, iwork, xwork);
  }
}

//
// Arguments    : const ::coder::array<double, 2U> &x
//                ::coder::array<int, 2U> &idx
// Return Type  : void
//
void sortIdx(const ::coder::array<double, 2U> &x, ::coder::array<int, 2U> &idx)
{
  array<int, 1U> iwork;
  int loop_ub;
  int n;
  int qEnd;
  n = x.size(1) + 1;
  idx.set_size(1, x.size(1));
  loop_ub = x.size(1);
  if (static_cast<int>(x.size(1) < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      idx[i] = 0;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      idx[i] = 0;
    }
  }
  if (x.size(1) != 0) {
    double d;
    int b_i;
    int k;
    iwork.set_size(x.size(1));
    loop_ub = x.size(1) - 1;
    for (k = 1; k <= loop_ub; k += 2) {
      d = x[k];
      if ((x[k - 1] <= d) || std::isnan(d)) {
        idx[k - 1] = k;
        idx[k] = k + 1;
      } else {
        idx[k - 1] = k + 1;
        idx[k] = k;
      }
    }
    if ((x.size(1) & 1) != 0) {
      idx[x.size(1) - 1] = x.size(1);
    }
    b_i = 2;
    while (b_i < n - 1) {
      int i2;
      int j;
      i2 = b_i << 1;
      j = 1;
      for (int pEnd{b_i + 1}; pEnd < n; pEnd = qEnd + b_i) {
        int kEnd;
        int p;
        int q;
        p = j;
        q = pEnd - 1;
        qEnd = j + i2;
        if (qEnd > n) {
          qEnd = n;
        }
        k = 0;
        kEnd = qEnd - j;
        while (k + 1 <= kEnd) {
          d = x[idx[q] - 1];
          loop_ub = idx[p - 1];
          if ((x[loop_ub - 1] <= d) || std::isnan(d)) {
            iwork[k] = loop_ub;
            p++;
            if (p == pEnd) {
              while (q + 1 < qEnd) {
                k++;
                iwork[k] = idx[q];
                q++;
              }
            }
          } else {
            iwork[k] = idx[q];
            q++;
            if (q + 1 == qEnd) {
              while (p < pEnd) {
                k++;
                iwork[k] = idx[p - 1];
                p++;
              }
            }
          }
          k++;
        }
        for (k = 0; k < kEnd; k++) {
          idx[(j + k) - 1] = iwork[k];
        }
        j = qEnd;
      }
      b_i = i2;
    }
  }
}

} // namespace internal
} // namespace coder

//
// File trailer for sortIdx.cpp
//
// [EOF]
//
