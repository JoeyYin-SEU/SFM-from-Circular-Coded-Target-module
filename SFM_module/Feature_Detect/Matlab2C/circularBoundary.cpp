//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: circularBoundary.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "circularBoundary.h"
#include "get_chessborad_pixel_rtwutil.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>

// Function Definitions
//
// Arguments    : const ::coder::array<float, 2U> &imagePoints
//                const ::coder::array<float, 2U> &b_I
//                ::coder::array<float, 2U> &points
// Return Type  : void
//
namespace coder {
namespace vision {
namespace internal {
namespace calibration {
namespace checkerboard {
void circularBoundary(const ::coder::array<float, 2U> &imagePoints,
                      const ::coder::array<float, 2U> &b_I,
                      ::coder::array<float, 2U> &points)
{
  static const float fv[36]{
      0.0F,           0.892784476F, 1.75687408F,  2.56449628F,  3.28969359F,
      3.90915751F,    4.40297747F,  4.7552824F,   4.95474863F,  4.99496555F,
      4.87463951F,    4.59763908F,  4.17286634F,  3.61397433F,  2.93892622F,
      2.16941881F,    1.33018422F,  0.44819653F,  -0.44819653F, -1.33018422F,
      -2.16941881F,   -2.93892622F, -3.61397433F, -4.17286634F, -4.59763908F,
      -4.87463951F,   -4.99496555F, -4.95474863F, -4.7552824F,  -4.40297747F,
      -3.90915751F,   -3.28969359F, -2.56449628F, -1.75687408F, -0.892784476F,
      -1.2246468E-15F};
  static const float fv1[36]{
      5.0F,         4.91964817F,   4.68117428F,  4.29224396F,  3.76535726F,
      3.11744905F,  2.36934328F,   1.54508495F,  0.671166301F, -0.224324152F,
      -1.11260462F, -1.9651252F,   -2.75448489F, -3.45531321F, -4.04508495F,
      -4.50484419F, -4.81981421F,  -4.97987127F, -4.97987127F, -4.81981421F,
      -4.50484419F, -4.04508495F,  -3.45531321F, -2.75448489F, -1.9651252F,
      -1.11260462F, -0.224324152F, 0.671166301F, 1.54508495F,  2.36934328F,
      3.11744905F,  3.76535726F,   4.29224396F,  4.68117428F,  4.91964817F,
      5.0F};
  array<float, 2U> pointDiff;
  array<float, 2U> validPoints;
  array<float, 1U> dists;
  array<int, 1U> r;
  array<int, 1U> r1;
  array<boolean_T, 1U> isValidPoint;
  float yBoundary[36];
  float b_imagePoints;
  float bsum;
  float f;
  float tmp2;
  float work;
  int idx[36];
  unsigned int imageSize[2];
  int siz[2];
  int acoef;
  int i;
  int i1;
  int k;
  int lastBlockLength;
  int n;
  int nblocks;
  unsigned int u2;
  unsigned int u3;
  boolean_T b;
  isValidPoint.set_size(imagePoints.size(0));
  acoef = imagePoints.size(0);
  if (static_cast<int>(imagePoints.size(0) < 3200)) {
    for (i = 0; i < acoef; i++) {
      isValidPoint[i] = false;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (i = 0; i < acoef; i++) {
      isValidPoint[i] = false;
    }
  }
  imageSize[0] = static_cast<unsigned int>(b_I.size(0));
  imageSize[1] = static_cast<unsigned int>(b_I.size(1));
  i1 = imagePoints.size(0);
  if (imagePoints.size(0) - 1 >= 0) {
    siz[0] = b_I.size(0);
  }
  if (static_cast<int>(imagePoints.size(0) * 36 < 3200)) {
    for (int b_i{0}; b_i < i1; b_i++) {
      unsigned int u;
      unsigned int u1;
      b_imagePoints = imagePoints[b_i];
      work = imagePoints[b_i + imagePoints.size(0)];
      u = imageSize[1];
      u1 = imageSize[0];
      for (k = 0; k < 36; k++) {
        idx[k] = static_cast<int>(std::fmax(std::fmin(std::round(work - fv[k]),
                                                      static_cast<float>(u1)),
                                            1.0F)) +
                 siz[0] * (static_cast<int>(std::fmax(
                               std::fmin(std::round(b_imagePoints - fv1[k]),
                                         static_cast<float>(u)),
                               1.0F)) -
                           1);
      }
      if (std::isnan(b_I[idx[0] - 1])) {
        tmp2 = 0.0F;
        n = 0;
      } else {
        tmp2 = b_I[idx[0] - 1];
        n = 1;
      }
      for (k = 0; k < 35; k++) {
        nblocks = idx[k + 1] - 1;
        if (!std::isnan(b_I[nblocks])) {
          tmp2 += b_I[nblocks];
          n++;
        }
      }
      tmp2 /= static_cast<float>(n);
      for (k = 0; k < 36; k++) {
        f = b_I[idx[k] - 1] - tmp2;
        if (std::isnan(f)) {
          f = rtNaNF;
        } else if (f < 0.0F) {
          f = -1.0F;
        } else {
          f = (f > 0.0F);
        }
        yBoundary[k] = f;
      }
      work = yBoundary[0];
      n = 0;
      for (k = 0; k < 35; k++) {
        tmp2 = work;
        bsum = yBoundary[k + 1];
        work = bsum;
        if (std::abs(bsum - tmp2) != 0.0F) {
          n++;
        }
      }
      isValidPoint[b_i] = (n > 2);
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads()                         \
                               : 32) private(f, tmp2, k, b, n, work, idx, u2,  \
                                             u3, yBoundary, b_imagePoints, i)

    for (int b_i = 0; b_i < i1; b_i++) {
      b_imagePoints = imagePoints[b_i];
      work = imagePoints[b_i + imagePoints.size(0)];
      u3 = imageSize[1];
      u2 = imageSize[0];
      n = siz[0];
      for (k = 0; k < 36; k++) {
        tmp2 = std::fmax(std::fmin(std::round(b_imagePoints - fv1[k]),
                                   static_cast<float>(u3)),
                         1.0F);
        f = std::fmax(
            std::fmin(std::round(work - fv[k]), static_cast<float>(u2)), 1.0F);
        idx[k] = static_cast<int>(f) + n * (static_cast<int>(tmp2) - 1);
      }
      if (std::isnan(b_I[idx[0] - 1])) {
        tmp2 = 0.0F;
        n = 0;
      } else {
        tmp2 = b_I[idx[0] - 1];
        n = 1;
      }
      for (k = 0; k < 35; k++) {
        i = idx[k + 1] - 1;
        b = std::isnan(b_I[i]);
        if (!b) {
          tmp2 += b_I[i];
          n++;
        }
      }
      tmp2 /= static_cast<float>(n);
      for (k = 0; k < 36; k++) {
        f = b_I[idx[k] - 1] - tmp2;
        if (std::isnan(f)) {
          f = rtNaNF;
        } else if (f < 0.0F) {
          f = -1.0F;
        } else {
          f = (f > 0.0F);
        }
        yBoundary[k] = f;
      }
      work = yBoundary[0];
      n = 0;
      for (k = 0; k < 35; k++) {
        tmp2 = work;
        work = yBoundary[k + 1];
        f = work - tmp2;
        f = std::abs(f);
        if (f != 0.0F) {
          n++;
        }
      }
      isValidPoint[b_i] = (n > 2);
    }
  }
  nblocks = isValidPoint.size(0) - 1;
  acoef = 0;
  for (lastBlockLength = 0; lastBlockLength <= nblocks; lastBlockLength++) {
    if (isValidPoint[lastBlockLength]) {
      acoef++;
    }
  }
  r.set_size(acoef);
  acoef = 0;
  for (lastBlockLength = 0; lastBlockLength <= nblocks; lastBlockLength++) {
    if (isValidPoint[lastBlockLength]) {
      r[acoef] = lastBlockLength + 1;
      acoef++;
    }
  }
  validPoints.set_size(r.size(0), 2);
  acoef = r.size(0);
  if (static_cast<int>((r.size(0) << 1) < 3200)) {
    for (i = 0; i < 2; i++) {
      for (n = 0; n < acoef; n++) {
        validPoints[n + validPoints.size(0) * i] =
            imagePoints[(r[n] + imagePoints.size(0) * i) - 1];
      }
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32) private(n)

    for (i = 0; i < 2; i++) {
      for (n = 0; n < acoef; n++) {
        validPoints[n + validPoints.size(0) * i] =
            imagePoints[(r[n] + imagePoints.size(0) * i) - 1];
      }
    }
  }
  points.set_size(0, 2);
  while (validPoints.size(0) != 0) {
    float x[2];
    pointDiff.set_size(validPoints.size(0), 2);
    acoef = (validPoints.size(0) != 1);
    for (int b_k{0}; b_k < 2; b_k++) {
      i1 = pointDiff.size(0) - 1;
      for (nblocks = 0; nblocks <= i1; nblocks++) {
        pointDiff[nblocks + pointDiff.size(0) * b_k] =
            validPoints[acoef * nblocks + validPoints.size(0) * b_k] -
            validPoints[validPoints.size(0) * b_k];
      }
    }
    dists.set_size(pointDiff.size(0));
    acoef = pointDiff.size(0);
    for (int b_k{0}; b_k < acoef; b_k++) {
      dists[b_k] =
          rt_hypotf_snf(pointDiff[b_k], pointDiff[b_k + pointDiff.size(0)]);
    }
    nblocks = dists.size(0) - 1;
    acoef = 0;
    for (lastBlockLength = 0; lastBlockLength <= nblocks; lastBlockLength++) {
      if (dists[lastBlockLength] < 5.0F) {
        acoef++;
      }
    }
    r1.set_size(acoef);
    acoef = 0;
    for (lastBlockLength = 0; lastBlockLength <= nblocks; lastBlockLength++) {
      if (dists[lastBlockLength] < 5.0F) {
        r1[acoef] = lastBlockLength + 1;
        acoef++;
      }
    }
    if (r1.size(0) == 0) {
      x[0] = 0.0F;
      x[1] = 0.0F;
    } else {
      if (r1.size(0) <= 1024) {
        acoef = r1.size(0);
        lastBlockLength = 0;
        nblocks = 1;
      } else {
        acoef = 1024;
        nblocks = static_cast<int>(static_cast<unsigned int>(r1.size(0)) >> 10);
        lastBlockLength = r1.size(0) - (nblocks << 10);
        if (lastBlockLength > 0) {
          nblocks++;
        } else {
          lastBlockLength = 1024;
        }
      }
      for (int xi{0}; xi < 2; xi++) {
        int xpageoffset;
        xpageoffset = xi * r1.size(0);
        x[xi] = validPoints[(r1[xpageoffset % r1.size(0)] +
                             validPoints.size(0) * (xpageoffset / r1.size(0))) -
                            1];
        for (int b_k{2}; b_k <= acoef; b_k++) {
          i1 = (xpageoffset + b_k) - 1;
          x[xi] += validPoints[(r1[i1 % r1.size(0)] +
                                validPoints.size(0) * (i1 / r1.size(0))) -
                               1];
        }
        for (int ib{2}; ib <= nblocks; ib++) {
          int hi;
          int xblockoffset;
          xblockoffset = xpageoffset + ((ib - 1) << 10);
          bsum =
              validPoints[(r1[xblockoffset % r1.size(0)] +
                           validPoints.size(0) * (xblockoffset / r1.size(0))) -
                          1];
          if (ib == nblocks) {
            hi = lastBlockLength;
          } else {
            hi = 1024;
          }
          for (int b_k{2}; b_k <= hi; b_k++) {
            i1 = (xblockoffset + b_k) - 1;
            bsum += validPoints[(r1[i1 % r1.size(0)] +
                                 validPoints.size(0) * (i1 / r1.size(0))) -
                                1];
          }
          x[xi] += bsum;
        }
      }
    }
    pointDiff.set_size(points.size(0) + 1, 2);
    acoef = points.size(0);
    for (int b_k{0}; b_k < 2; b_k++) {
      x[b_k] = std::round(x[b_k] / static_cast<float>(r1.size(0)));
      for (i1 = 0; i1 < acoef; i1++) {
        pointDiff[i1 + pointDiff.size(0) * b_k] =
            points[i1 + points.size(0) * b_k];
      }
    }
    pointDiff[points.size(0)] = x[0];
    pointDiff[points.size(0) + pointDiff.size(0)] = x[1];
    points.set_size(pointDiff.size(0), 2);
    acoef = pointDiff.size(0) * 2;
    for (i1 = 0; i1 < acoef; i1++) {
      points[i1] = pointDiff[i1];
    }
    isValidPoint.set_size(dists.size(0));
    acoef = dists.size(0);
    for (i1 = 0; i1 < acoef; i1++) {
      isValidPoint[i1] = (dists[i1] < 5.0F);
    }
    acoef = validPoints.size(0);
    nblocks = 0;
    i1 = isValidPoint.size(0);
    for (int b_k{0}; b_k < i1; b_k++) {
      nblocks += isValidPoint[b_k];
    }
    nblocks = validPoints.size(0) - nblocks;
    lastBlockLength = 0;
    for (int b_k{0}; b_k < acoef; b_k++) {
      if ((b_k + 1 > isValidPoint.size(0)) || (!isValidPoint[b_k])) {
        validPoints[lastBlockLength] = validPoints[b_k];
        validPoints[lastBlockLength + validPoints.size(0)] =
            validPoints[b_k + validPoints.size(0)];
        lastBlockLength++;
      }
    }
    if (nblocks < 1) {
      acoef = 0;
    } else {
      acoef = nblocks;
    }
    for (i1 = 0; i1 < 2; i1++) {
      for (nblocks = 0; nblocks < acoef; nblocks++) {
        validPoints[nblocks + acoef * i1] =
            validPoints[nblocks + validPoints.size(0) * i1];
      }
    }
    validPoints.set_size(acoef, 2);
  }
}

} // namespace checkerboard
} // namespace calibration
} // namespace internal
} // namespace vision
} // namespace coder

//
// File trailer for circularBoundary.cpp
//
// [EOF]
//
