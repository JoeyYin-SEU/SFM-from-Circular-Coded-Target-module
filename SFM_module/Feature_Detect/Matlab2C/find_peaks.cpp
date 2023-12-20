//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: find_peaks.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "find_peaks.h"
#include "algbwmorph.h"
#include "get_chessborad_pixel_rtwutil.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "libmwimregionalmax.h"
#include "omp.h"
#include <cmath>

// Function Definitions
//
// Arguments    : const ::coder::array<float, 2U> &metric
//                double quality
//                ::coder::array<float, 2U> &loc
// Return Type  : void
//
namespace coder {
namespace vision {
namespace internal {
namespace calibration {
namespace checkerboard {
void find_peaks(const ::coder::array<float, 2U> &metric, double quality,
                ::coder::array<float, 2U> &loc)
{
  array<int, 1U> ii;
  array<int, 1U> vk;
  array<boolean_T, 2U> bw;
  array<boolean_T, 2U> last_aout;
  float maxMetric;
  int idx;
  int k;
  int last;
  boolean_T exitg1;
  last = metric.size(0) * metric.size(1);
  if (metric.size(0) * metric.size(1) <= 2) {
    if (metric.size(0) * metric.size(1) == 1) {
      maxMetric = metric[0];
    } else if ((metric[0] < metric[metric.size(0) * metric.size(1) - 1]) ||
               (std::isnan(metric[0]) &&
                (!std::isnan(metric[metric.size(0) * metric.size(1) - 1])))) {
      maxMetric = metric[metric.size(0) * metric.size(1) - 1];
    } else {
      maxMetric = metric[0];
    }
  } else {
    if (!std::isnan(metric[0])) {
      idx = 1;
    } else {
      idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= last)) {
        if (!std::isnan(metric[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }
    if (idx == 0) {
      maxMetric = metric[0];
    } else {
      maxMetric = metric[idx - 1];
      idx++;
      for (k = idx; k <= last; k++) {
        if (maxMetric < metric[k - 1]) {
          maxMetric = metric[k - 1];
        }
      }
    }
  }
  if (maxMetric <= 4.94065645841247E-324) {
    loc.set_size(0, 2);
  } else {
    double connSizeT[2];
    double imSizeT[2];
    boolean_T conn[9];
    bw.set_size(metric.size(0), metric.size(1));
    imSizeT[0] = metric.size(0);
    imSizeT[1] = metric.size(1);
    for (idx = 0; idx < 9; idx++) {
      conn[idx] = true;
    }
    connSizeT[0] = 3.0;
    connSizeT[1] = 3.0;
    imregionalmax_real32(&metric[0], &bw[0], 2.0, &imSizeT[0], &conn[0], 2.0,
                         &connSizeT[0]);
    maxMetric *= static_cast<float>(quality);
    k = metric.size(0) * metric.size(1) - 1;
    last = 0;
    for (idx = 0; idx <= k; idx++) {
      if (metric[idx] < maxMetric) {
        last++;
      }
    }
    ii.set_size(last);
    last = 0;
    for (idx = 0; idx <= k; idx++) {
      if (metric[idx] < maxMetric) {
        ii[last] = idx + 1;
        last++;
      }
    }
    k = ii.size(0);
    for (idx = 0; idx < k; idx++) {
      bw[ii[idx] - 1] = false;
    }
    if ((bw.size(0) != 0) && (bw.size(1) != 0)) {
      boolean_T p;
      do {
        last_aout.set_size(bw.size(0), bw.size(1));
        k = bw.size(0) * bw.size(1);
        for (idx = 0; idx < k; idx++) {
          last_aout[idx] = bw[idx];
        }
        images::internal::bwmorphApplyOnce(bw);
        p = false;
        if ((last_aout.size(0) == bw.size(0)) &&
            (last_aout.size(1) == bw.size(1))) {
          p = true;
        }
        if (p && ((last_aout.size(0) != 0) && (last_aout.size(1) != 0)) &&
            ((bw.size(0) != 0) && (bw.size(1) != 0))) {
          k = 0;
          exitg1 = false;
          while ((!exitg1) && (k <= bw.size(0) * bw.size(1) - 1)) {
            if (last_aout[k] != bw[k]) {
              p = false;
              exitg1 = true;
            } else {
              k++;
            }
          }
        }
      } while (!p);
      //  the output is not changing anymore
    }
    k = bw.size(1);
    if (static_cast<int>(bw.size(1) < 3200)) {
      for (int i{0}; i < k; i++) {
        bw[bw.size(0) * i] = false;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < k; i++) {
        bw[bw.size(0) * i] = false;
      }
    }
    last = bw.size(0);
    k = bw.size(1);
    if (static_cast<int>(bw.size(1) < 3200)) {
      for (int i{0}; i < k; i++) {
        bw[(last + bw.size(0) * i) - 1] = false;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < k; i++) {
        bw[(last + bw.size(0) * i) - 1] = false;
      }
    }
    k = bw.size(0);
    if (static_cast<int>(bw.size(0) < 3200)) {
      for (int i{0}; i < k; i++) {
        bw[i] = false;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < k; i++) {
        bw[i] = false;
      }
    }
    last = bw.size(1);
    k = bw.size(0);
    if (static_cast<int>(bw.size(0) < 3200)) {
      for (int i{0}; i < k; i++) {
        bw[i + bw.size(0) * (last - 1)] = false;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < k; i++) {
        bw[i + bw.size(0) * (last - 1)] = false;
      }
    }
    last = bw.size(0) * bw.size(1);
    idx = 0;
    ii.set_size(last);
    k = 0;
    exitg1 = false;
    while ((!exitg1) && (k <= last - 1)) {
      if (bw[k]) {
        idx++;
        ii[idx - 1] = k + 1;
        if (idx >= last) {
          exitg1 = true;
        } else {
          k++;
        }
      } else {
        k++;
      }
    }
    if (last == 1) {
      if (idx == 0) {
        ii.set_size(0);
      }
    } else {
      if (idx < 1) {
        idx = 0;
      }
      ii.set_size(idx);
    }
    loc.set_size(ii.size(0), 2);
    k = ii.size(0) << 1;
    if (static_cast<int>(k < 3200)) {
      for (int i{0}; i < k; i++) {
        loc[i] = 0.0F;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < k; i++) {
        loc[i] = 0.0F;
      }
    }
    connSizeT[0] = metric.size(0);
    k = ii.size(0);
    if (static_cast<int>(ii.size(0) < 3200)) {
      for (int i{0}; i < k; i++) {
        ii[i] = ii[i] - 1;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < k; i++) {
        ii[i] = ii[i] - 1;
      }
    }
    last = static_cast<int>(connSizeT[0]);
    vk.set_size(ii.size(0));
    k = ii.size(0);
    if (static_cast<int>(ii.size(0) < 3200)) {
      for (int i{0}; i < k; i++) {
        vk[i] = div_s32(ii[i], last);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < k; i++) {
        vk[i] = div_s32(ii[i], last);
      }
    }
    last = static_cast<int>(connSizeT[0]);
    k = ii.size(0);
    if (static_cast<int>(ii.size(0) < 3200)) {
      for (int i{0}; i < k; i++) {
        ii[i] = ii[i] - vk[i] * last;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < k; i++) {
        ii[i] = ii[i] - vk[i] * last;
      }
    }
    k = ii.size(0);
    if (static_cast<int>(ii.size(0) < 3200)) {
      for (int i{0}; i < k; i++) {
        loc[i + loc.size(0)] = static_cast<float>(ii[i] + 1);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < k; i++) {
        loc[i + loc.size(0)] = static_cast<float>(ii[i] + 1);
      }
    }
    k = vk.size(0);
    if (static_cast<int>(vk.size(0) < 3200)) {
      for (int i{0}; i < k; i++) {
        loc[i] = static_cast<float>(vk[i] + 1);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < k; i++) {
        loc[i] = static_cast<float>(vk[i] + 1);
      }
    }
  }
}

} // namespace checkerboard
} // namespace calibration
} // namespace internal
} // namespace vision
} // namespace coder

//
// File trailer for find_peaks.cpp
//
// [EOF]
//
