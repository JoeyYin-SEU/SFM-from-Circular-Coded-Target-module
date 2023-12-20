//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: ismember.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "ismember.h"
#include "rt_nonfinite.h"
#include "sort.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>

// Function Declarations
namespace coder {
static int bsearchni(int k, const ::coder::array<double, 2U> &x,
                     const ::coder::array<double, 2U> &s);

static int bsearchni(int k, const ::coder::array<double, 2U> &x,
                     const ::coder::array<double, 1U> &s);

} // namespace coder

// Function Definitions
//
// Arguments    : int k
//                const ::coder::array<double, 2U> &x
//                const ::coder::array<double, 2U> &s
// Return Type  : int
//
namespace coder {
static int bsearchni(int k, const ::coder::array<double, 2U> &x,
                     const ::coder::array<double, 2U> &s)
{
  double b_x;
  int idx;
  int ihi;
  int ilo;
  boolean_T exitg1;
  b_x = x[k - 1];
  ihi = s.size(1);
  idx = 0;
  ilo = 1;
  exitg1 = false;
  while ((!exitg1) && (ihi >= ilo)) {
    int imid;
    imid = ((ilo >> 1) + (ihi >> 1)) - 1;
    if (((ilo & 1) == 1) && ((ihi & 1) == 1)) {
      imid++;
    }
    if (b_x == s[imid]) {
      idx = imid + 1;
      exitg1 = true;
    } else {
      boolean_T p;
      if (std::isnan(s[imid])) {
        p = !std::isnan(b_x);
      } else if (std::isnan(b_x)) {
        p = false;
      } else {
        p = (b_x < s[imid]);
      }
      if (p) {
        ihi = imid;
      } else {
        ilo = imid + 2;
      }
    }
  }
  if (idx > 0) {
    idx--;
    while ((idx > 0) && (b_x == s[idx - 1])) {
      idx--;
    }
    idx++;
  }
  return idx;
}

//
// Arguments    : int k
//                const ::coder::array<double, 2U> &x
//                const ::coder::array<double, 1U> &s
// Return Type  : int
//
static int bsearchni(int k, const ::coder::array<double, 2U> &x,
                     const ::coder::array<double, 1U> &s)
{
  double b_x;
  int idx;
  int ihi;
  int ilo;
  boolean_T exitg1;
  b_x = x[k - 1];
  ihi = s.size(0);
  idx = 0;
  ilo = 1;
  exitg1 = false;
  while ((!exitg1) && (ihi >= ilo)) {
    int imid;
    imid = ((ilo >> 1) + (ihi >> 1)) - 1;
    if (((ilo & 1) == 1) && ((ihi & 1) == 1)) {
      imid++;
    }
    if (b_x == s[imid]) {
      idx = imid + 1;
      exitg1 = true;
    } else {
      boolean_T p;
      if (std::isnan(s[imid])) {
        p = !std::isnan(b_x);
      } else if (std::isnan(b_x)) {
        p = false;
      } else {
        p = (b_x < s[imid]);
      }
      if (p) {
        ihi = imid;
      } else {
        ilo = imid + 2;
      }
    }
  }
  if (idx > 0) {
    idx--;
    while ((idx > 0) && (b_x == s[idx - 1])) {
      idx--;
    }
    idx++;
  }
  return idx;
}

//
// Arguments    : const ::coder::array<double, 2U> &a
//                const ::coder::array<double, 2U> &s
//                ::coder::array<boolean_T, 2U> &tf
// Return Type  : void
//
void isMember(const ::coder::array<double, 2U> &a,
              const ::coder::array<double, 2U> &s,
              ::coder::array<boolean_T, 2U> &tf)
{
  array<double, 1U> ss;
  array<int, 1U> b_ss;
  int k;
  int n;
  int na;
  int ns;
  int pmax;
  int pmin;
  boolean_T exitg1;
  boolean_T guard1{false};
  na = a.size(1);
  ns = s.size(1);
  pmax = a.size(1);
  tf.set_size(1, a.size(1));
  pmin = a.size(1);
  if (static_cast<int>(a.size(1) < 3200)) {
    for (n = 0; n < pmax; n++) {
      tf[n] = false;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (n = 0; n < pmin; n++) {
      tf[n] = false;
    }
  }
  guard1 = false;
  if (s.size(1) <= 4) {
    guard1 = true;
  } else {
    pmax = 31;
    pmin = 0;
    exitg1 = false;
    while ((!exitg1) && (pmax - pmin > 1)) {
      int p;
      int pow2p;
      p = (pmin + pmax) >> 1;
      pow2p = 1 << p;
      if (pow2p == ns) {
        pmax = p;
        exitg1 = true;
      } else if (pow2p > ns) {
        pmax = p;
      } else {
        pmin = p;
      }
    }
    if (a.size(1) <= pmax + 4) {
      guard1 = true;
    } else {
      boolean_T y;
      y = true;
      pmax = 0;
      exitg1 = false;
      while ((!exitg1) && (pmax <= s.size(1) - 2)) {
        double v_idx_1;
        v_idx_1 = s[pmax + 1];
        if ((s[pmax] <= v_idx_1) || std::isnan(v_idx_1)) {
          pmax++;
        } else {
          y = false;
          exitg1 = true;
        }
      }
      if (!y) {
        ss.set_size(s.size(1));
        pmax = s.size(1);
        if (static_cast<int>(s.size(1) < 3200)) {
          for (n = 0; n < pmax; n++) {
            ss[n] = s[n];
          }
        } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

          for (n = 0; n < pmax; n++) {
            ss[n] = s[n];
          }
        }
        internal::sort(ss, b_ss);
        pmax = a.size(1) - 1;
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32) private(n)

        for (k = 0; k <= pmax; k++) {
          n = bsearchni(k + 1, a, ss);
          if (n > 0) {
            tf[k] = true;
          }
        }
      } else {
        pmax = a.size(1) - 1;
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32) private(n)

        for (k = 0; k <= pmax; k++) {
          n = bsearchni(k + 1, a, s);
          if (n > 0) {
            tf[k] = true;
          }
        }
      }
    }
  }
  if (guard1) {
    if (static_cast<int>(na * 10 < 3200)) {
      for (n = 0; n < na; n++) {
        k = 0;
        exitg1 = false;
        while ((!exitg1) && (k <= ns - 1)) {
          if (a[n] == s[k]) {
            tf[n] = true;
            exitg1 = true;
          } else {
            k++;
          }
        }
      }
    } else {
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(k, exitg1)

      for (n = 0; n < na; n++) {
        k = 0;
        exitg1 = false;
        while ((!exitg1) && (k <= ns - 1)) {
          if (a[n] == s[k]) {
            tf[n] = true;
            exitg1 = true;
          } else {
            k++;
          }
        }
      }
    }
  }
}

} // namespace coder

//
// File trailer for ismember.cpp
//
// [EOF]
//
