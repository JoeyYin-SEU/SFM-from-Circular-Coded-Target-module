//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: combineVectorElements.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "combineVectorElements.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"

// Function Definitions
//
// Arguments    : const ::coder::array<boolean_T, 2U> &x
//                ::coder::array<int, 2U> &y
// Return Type  : void
//
namespace coder {
void b_combineVectorElements(const ::coder::array<boolean_T, 2U> &x,
                             ::coder::array<int, 2U> &y)
{
  int k;
  int vlen;
  int xpageoffset;
  vlen = x.size(0);
  if ((x.size(0) == 0) || (x.size(1) == 0)) {
    y.set_size(1, x.size(1));
    vlen = x.size(1);
    if (static_cast<int>(x.size(1) < 3200)) {
      for (xpageoffset = 0; xpageoffset < vlen; xpageoffset++) {
        y[xpageoffset] = 0;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (xpageoffset = 0; xpageoffset < vlen; xpageoffset++) {
        y[xpageoffset] = 0;
      }
    }
  } else {
    int npages;
    npages = x.size(1);
    y.set_size(1, x.size(1));
    if (static_cast<int>(x.size(1) * (x.size(0) - 1) < 3200)) {
      for (int i{0}; i < npages; i++) {
        xpageoffset = i * x.size(0);
        y[i] = x[xpageoffset];
        for (k = 2; k <= vlen; k++) {
          y[i] = y[i] + x[(xpageoffset + k) - 1];
        }
      }
    } else {
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(k, xpageoffset)

      for (int i = 0; i < npages; i++) {
        xpageoffset = i * x.size(0);
        y[i] = x[xpageoffset];
        for (k = 2; k <= vlen; k++) {
          y[i] = y[i] + x[(xpageoffset + k) - 1];
        }
      }
    }
  }
}

//
// Arguments    : const ::coder::array<double, 1U> &x
// Return Type  : double
//
double combineVectorElements(const ::coder::array<double, 1U> &x)
{
  double y;
  if (x.size(0) == 0) {
    y = 0.0;
  } else {
    int firstBlockLength;
    int lastBlockLength;
    int nblocks;
    if (x.size(0) <= 1024) {
      firstBlockLength = x.size(0);
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      firstBlockLength = 1024;
      nblocks = static_cast<int>(static_cast<unsigned int>(x.size(0)) >> 10);
      lastBlockLength = x.size(0) - (nblocks << 10);
      if (lastBlockLength > 0) {
        nblocks++;
      } else {
        lastBlockLength = 1024;
      }
    }
    y = x[0];
    for (int k{2}; k <= firstBlockLength; k++) {
      y += x[k - 1];
    }
    for (int ib{2}; ib <= nblocks; ib++) {
      double bsum;
      int hi;
      firstBlockLength = (ib - 1) << 10;
      bsum = x[firstBlockLength];
      if (ib == nblocks) {
        hi = lastBlockLength;
      } else {
        hi = 1024;
      }
      for (int k{2}; k <= hi; k++) {
        bsum += x[(firstBlockLength + k) - 1];
      }
      y += bsum;
    }
  }
  return y;
}

//
// Arguments    : const ::coder::array<double, 2U> &x
// Return Type  : double
//
double combineVectorElements(const ::coder::array<double, 2U> &x)
{
  double y;
  if (x.size(1) == 0) {
    y = 0.0;
  } else {
    int firstBlockLength;
    int lastBlockLength;
    int nblocks;
    if (x.size(1) <= 1024) {
      firstBlockLength = x.size(1);
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      firstBlockLength = 1024;
      nblocks = static_cast<int>(static_cast<unsigned int>(x.size(1)) >> 10);
      lastBlockLength = x.size(1) - (nblocks << 10);
      if (lastBlockLength > 0) {
        nblocks++;
      } else {
        lastBlockLength = 1024;
      }
    }
    y = x[0];
    for (int k{2}; k <= firstBlockLength; k++) {
      y += x[k - 1];
    }
    for (int ib{2}; ib <= nblocks; ib++) {
      double bsum;
      int hi;
      firstBlockLength = (ib - 1) << 10;
      bsum = x[firstBlockLength];
      if (ib == nblocks) {
        hi = lastBlockLength;
      } else {
        hi = 1024;
      }
      for (int k{2}; k <= hi; k++) {
        bsum += x[(firstBlockLength + k) - 1];
      }
      y += bsum;
    }
  }
  return y;
}

} // namespace coder

//
// File trailer for combineVectorElements.cpp
//
// [EOF]
//
