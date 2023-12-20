//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: minOrMax.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "minOrMax.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>

// Function Definitions
//
// Arguments    : const float x[3]
// Return Type  : float
//
namespace coder {
namespace internal {
float maximum(const float x[3])
{
  float ex;
  int idx;
  int k;
  if (!std::isnan(x[0])) {
    idx = 1;
  } else {
    boolean_T exitg1;
    idx = 0;
    k = 2;
    exitg1 = false;
    while ((!exitg1) && (k <= 3)) {
      if (!std::isnan(x[k - 1])) {
        idx = k;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }
  if (idx == 0) {
    ex = x[0];
  } else {
    ex = x[idx - 1];
    idx++;
    for (k = idx; k < 4; k++) {
      float f;
      f = x[k - 1];
      if (ex < f) {
        ex = f;
      }
    }
  }
  return ex;
}

//
// Arguments    : const ::coder::array<double, 1U> &x
// Return Type  : double
//
double maximum(const ::coder::array<double, 1U> &x)
{
  double ex;
  int last;
  last = x.size(0);
  if (x.size(0) <= 2) {
    if (x.size(0) == 1) {
      ex = x[0];
    } else if ((x[0] < x[x.size(0) - 1]) ||
               (std::isnan(x[0]) && (!std::isnan(x[x.size(0) - 1])))) {
      ex = x[x.size(0) - 1];
    } else {
      ex = x[0];
    }
  } else {
    int idx;
    int k;
    if (!std::isnan(x[0])) {
      idx = 1;
    } else {
      boolean_T exitg1;
      idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= last)) {
        if (!std::isnan(x[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }
    if (idx == 0) {
      ex = x[0];
    } else {
      ex = x[idx - 1];
      idx++;
      for (k = idx; k <= last; k++) {
        double d;
        d = x[k - 1];
        if (ex < d) {
          ex = d;
        }
      }
    }
  }
  return ex;
}

//
// Arguments    : const ::coder::array<double, 2U> &x
//                ::coder::array<double, 2U> &ex
// Return Type  : void
//
void maximum(const ::coder::array<double, 2U> &x,
             ::coder::array<double, 2U> &ex)
{
  double a;
  double b;
  int i;
  int m;
  int n;
  boolean_T p;
  m = x.size(0);
  n = x.size(1);
  ex.set_size(1, x.size(1));
  if (x.size(1) >= 1) {
    if (static_cast<int>(x.size(1) * (x.size(0) - 1) < 3200)) {
      for (int j{0}; j < n; j++) {
        ex[j] = x[x.size(0) * j];
        for (i = 2; i <= m; i++) {
          a = ex[j];
          b = x[(i + x.size(0) * j) - 1];
          if (std::isnan(b)) {
            p = false;
          } else if (std::isnan(a)) {
            p = true;
          } else {
            p = (a < b);
          }
          if (p) {
            ex[j] = b;
          }
        }
      }
    } else {
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(p, b, a, i)

      for (int j = 0; j < n; j++) {
        ex[j] = x[x.size(0) * j];
        for (i = 2; i <= m; i++) {
          a = ex[j];
          b = x[(i + x.size(0) * j) - 1];
          if (std::isnan(b)) {
            p = false;
          } else if (std::isnan(a)) {
            p = true;
          } else {
            p = (a < b);
          }
          if (p) {
            ex[j] = b;
          }
        }
      }
    }
  }
}

//
// Arguments    : const ::coder::array<float, 1U> &x
//                float *ex
//                int *idx
// Return Type  : void
//
void minimum(const ::coder::array<float, 1U> &x, float *ex, int *idx)
{
  int last;
  last = x.size(0);
  if (x.size(0) <= 2) {
    if (x.size(0) == 1) {
      *ex = x[0];
      *idx = 1;
    } else if ((x[0] > x[x.size(0) - 1]) ||
               (std::isnan(x[0]) && (!std::isnan(x[x.size(0) - 1])))) {
      *ex = x[x.size(0) - 1];
      *idx = x.size(0);
    } else {
      *ex = x[0];
      *idx = 1;
    }
  } else {
    int k;
    if (!std::isnan(x[0])) {
      *idx = 1;
    } else {
      boolean_T exitg1;
      *idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= last)) {
        if (!std::isnan(x[k - 1])) {
          *idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }
    if (*idx == 0) {
      *ex = x[0];
      *idx = 1;
    } else {
      int i;
      *ex = x[*idx - 1];
      i = *idx + 1;
      for (k = i; k <= last; k++) {
        float f;
        f = x[k - 1];
        if (*ex > f) {
          *ex = f;
          *idx = k;
        }
      }
    }
  }
}

//
// Arguments    : const ::coder::array<float, 1U> &x
// Return Type  : float
//
float minimum(const ::coder::array<float, 1U> &x)
{
  float ex;
  int last;
  last = x.size(0);
  if (x.size(0) <= 2) {
    if (x.size(0) == 1) {
      ex = x[0];
    } else if ((x[0] > x[x.size(0) - 1]) ||
               (std::isnan(x[0]) && (!std::isnan(x[x.size(0) - 1])))) {
      ex = x[x.size(0) - 1];
    } else {
      ex = x[0];
    }
  } else {
    int idx;
    int k;
    if (!std::isnan(x[0])) {
      idx = 1;
    } else {
      boolean_T exitg1;
      idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= last)) {
        if (!std::isnan(x[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }
    if (idx == 0) {
      ex = x[0];
    } else {
      ex = x[idx - 1];
      idx++;
      for (k = idx; k <= last; k++) {
        float f;
        f = x[k - 1];
        if (ex > f) {
          ex = f;
        }
      }
    }
  }
  return ex;
}

} // namespace internal
} // namespace coder

//
// File trailer for minOrMax.cpp
//
// [EOF]
//
