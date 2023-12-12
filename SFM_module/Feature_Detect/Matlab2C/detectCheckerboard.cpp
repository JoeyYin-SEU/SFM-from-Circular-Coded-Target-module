//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: detectCheckerboard.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

// Include Files
#include "detectCheckerboard.h"
#include "Checkerboard.h"
#include "flip.h"
#include "isequal.h"
#include "mean.h"
#include "rot90.h"
#include "rt_nonfinite.h"
#include "sort.h"
#include "coder_array.h"
#include "omp.h"
#include "rt_defines.h"
#include <cmath>

// Function Declarations
namespace coder {
namespace vision {
namespace internal {
namespace calibration {
namespace checkerboard {
static boolean_T isUpperLeftBlack(const Checkerboard *b_this,
                                  const ::coder::array<float, 2U> &b_I);

static void poly2RectMask(double b_X[4], double Y[4], double height,
                          double width, ::coder::array<boolean_T, 2U> &mask);

} // namespace checkerboard
} // namespace calibration
} // namespace internal
} // namespace vision
} // namespace coder
static float rt_atan2f_snf(float u0, float u1);

// Function Definitions
//
// Arguments    : const Checkerboard *b_this
//                const ::coder::array<float, 2U> &b_I
// Return Type  : boolean_T
//
namespace coder {
namespace vision {
namespace internal {
namespace calibration {
namespace checkerboard {
static boolean_T isUpperLeftBlack(const Checkerboard *b_this,
                                  const ::coder::array<float, 2U> &b_I)
{
  array<float, 1U> c_I;
  array<float, 1U> d_I;
  array<int, 1U> r;
  array<int, 1U> r1;
  array<boolean_T, 2U> nextSquareMask;
  array<boolean_T, 2U> upperLeftMask;
  double c_this[4];
  double d_this[4];
  int end;
  int trueCount;
  boolean_T tf;
  c_this[0] = b_this->BoardCoords[0];
  c_this[1] = b_this->BoardCoords[b_this->BoardCoords.size(0)];
  c_this[2] = b_this->BoardCoords[b_this->BoardCoords.size(0) + 1];
  c_this[3] = b_this->BoardCoords[1];
  d_this[0] = b_this->BoardCoords[b_this->BoardCoords.size(0) *
                                  b_this->BoardCoords.size(1)];
  d_this[1] = b_this->BoardCoords[b_this->BoardCoords.size(0) +
                                  b_this->BoardCoords.size(0) *
                                      b_this->BoardCoords.size(1)];
  d_this[2] = b_this->BoardCoords[(b_this->BoardCoords.size(0) +
                                   b_this->BoardCoords.size(0) *
                                       b_this->BoardCoords.size(1)) +
                                  1];
  d_this[3] = b_this->BoardCoords[b_this->BoardCoords.size(0) *
                                      b_this->BoardCoords.size(1) +
                                  1];
  poly2RectMask(c_this, d_this, static_cast<double>(b_I.size(0)),
                static_cast<double>(b_I.size(1)), upperLeftMask);
  c_this[0] = b_this->BoardCoords[b_this->BoardCoords.size(0)];
  c_this[1] = b_this->BoardCoords[b_this->BoardCoords.size(0) * 2];
  c_this[2] = b_this->BoardCoords[b_this->BoardCoords.size(0) * 2 + 1];
  c_this[3] = b_this->BoardCoords[b_this->BoardCoords.size(0) + 1];
  d_this[0] = b_this->BoardCoords[b_this->BoardCoords.size(0) +
                                  b_this->BoardCoords.size(0) *
                                      b_this->BoardCoords.size(1)];
  d_this[1] = b_this->BoardCoords[b_this->BoardCoords.size(0) * 2 +
                                  b_this->BoardCoords.size(0) *
                                      b_this->BoardCoords.size(1)];
  d_this[2] = b_this->BoardCoords[(b_this->BoardCoords.size(0) * 2 +
                                   b_this->BoardCoords.size(0) *
                                       b_this->BoardCoords.size(1)) +
                                  1];
  d_this[3] = b_this->BoardCoords[(b_this->BoardCoords.size(0) +
                                   b_this->BoardCoords.size(0) *
                                       b_this->BoardCoords.size(1)) +
                                  1];
  poly2RectMask(c_this, d_this, static_cast<double>(b_I.size(0)),
                static_cast<double>(b_I.size(1)), nextSquareMask);
  end = upperLeftMask.size(0) * upperLeftMask.size(1) - 1;
  trueCount = 0;
  for (int i{0}; i <= end; i++) {
    if (upperLeftMask[i]) {
      trueCount++;
    }
  }
  r.set_size(trueCount);
  trueCount = 0;
  for (int i{0}; i <= end; i++) {
    if (upperLeftMask[i]) {
      r[trueCount] = i + 1;
      trueCount++;
    }
  }
  end = nextSquareMask.size(0) * nextSquareMask.size(1) - 1;
  trueCount = 0;
  for (int i{0}; i <= end; i++) {
    if (nextSquareMask[i]) {
      trueCount++;
    }
  }
  r1.set_size(trueCount);
  trueCount = 0;
  for (int i{0}; i <= end; i++) {
    if (nextSquareMask[i]) {
      r1[trueCount] = i + 1;
      trueCount++;
    }
  }
  c_I.set_size(r.size(0));
  trueCount = r.size(0);
  if (static_cast<int>(r.size(0) < 3200)) {
    for (int b_i{0}; b_i < trueCount; b_i++) {
      c_I[b_i] = b_I[r[b_i] - 1];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int b_i = 0; b_i < trueCount; b_i++) {
      c_I[b_i] = b_I[r[b_i] - 1];
    }
  }
  d_I.set_size(r1.size(0));
  trueCount = r1.size(0);
  if (static_cast<int>(r1.size(0) < 3200)) {
    for (int b_i{0}; b_i < trueCount; b_i++) {
      d_I[b_i] = b_I[r1[b_i] - 1];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int b_i = 0; b_i < trueCount; b_i++) {
      d_I[b_i] = b_I[r1[b_i] - 1];
    }
  }
  tf = (mean(c_I) < mean(d_I));
  return tf;
}

//
// Arguments    : double b_X[4]
//                double Y[4]
//                double height
//                double width
//                ::coder::array<boolean_T, 2U> &mask
// Return Type  : void
//
static void poly2RectMask(double b_X[4], double Y[4], double height,
                          double width, ::coder::array<boolean_T, 2U> &mask)
{
  int iv[2];
  int b_loop_ub;
  int i1;
  int i2;
  int i3;
  int i4;
  int loop_ub;
  ::coder::internal::sort(b_X);
  ::coder::internal::sort(Y);
  mask.set_size(static_cast<int>(height), static_cast<int>(width));
  loop_ub = static_cast<int>(height) * static_cast<int>(width);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      mask[i] = false;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      mask[i] = false;
    }
  }
  if (Y[1] > Y[2]) {
    i1 = 0;
    loop_ub = 0;
  } else {
    i1 = static_cast<int>(Y[1]) - 1;
    loop_ub = static_cast<int>(Y[2]);
  }
  if (b_X[1] > b_X[2]) {
    i2 = 0;
    i3 = 0;
  } else {
    i2 = static_cast<int>(b_X[1]) - 1;
    i3 = static_cast<int>(b_X[2]);
  }
  iv[0] = loop_ub - i1;
  iv[1] = i3 - i2;
  loop_ub = iv[1];
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(i4, b_loop_ub)

  for (int i = 0; i < loop_ub; i++) {
    b_loop_ub = iv[0];
    for (i4 = 0; i4 < b_loop_ub; i4++) {
      mask[(i1 + i4) + mask.size(0) * (i2 + i)] = true;
    }
  }
}

//
// Arguments    : float u0
//                float u1
// Return Type  : float
//
} // namespace checkerboard
} // namespace calibration
} // namespace internal
} // namespace vision
} // namespace coder
static float rt_atan2f_snf(float u0, float u1)
{
  float y;
  if (std::isnan(u0) || std::isnan(u1)) {
    y = rtNaNF;
  } else if (std::isinf(u0) && std::isinf(u1)) {
    int i;
    int i1;
    if (u0 > 0.0F) {
      i = 1;
    } else {
      i = -1;
    }
    if (u1 > 0.0F) {
      i1 = 1;
    } else {
      i1 = -1;
    }
    y = std::atan2(static_cast<float>(i), static_cast<float>(i1));
  } else if (u1 == 0.0F) {
    if (u0 > 0.0F) {
      y = RT_PIF / 2.0F;
    } else if (u0 < 0.0F) {
      y = -(RT_PIF / 2.0F);
    } else {
      y = 0.0F;
    }
  } else {
    y = std::atan2(u0, u1);
  }
  return y;
}

//
// Arguments    : const ::coder::array<float, 2U> &points
//                const ::coder::array<float, 1U> &scores
//                const ::coder::array<float, 2U> &Ix2
//                const ::coder::array<float, 2U> &Iy2
//                const ::coder::array<float, 2U> &Ixy
//                double theta
//                boolean_T highDistortion
//                boolean_T usePartial
//                Checkerboard *iobj_0
// Return Type  : Checkerboard *
//
namespace coder {
namespace vision {
namespace internal {
namespace calibration {
namespace checkerboard {
Checkerboard *growCheckerboard(const ::coder::array<float, 2U> &points,
                               const ::coder::array<float, 1U> &scores,
                               const ::coder::array<float, 2U> &Ix2,
                               const ::coder::array<float, 2U> &Iy2,
                               const ::coder::array<float, 2U> &Ixy,
                               double theta, boolean_T highDistortion,
                               boolean_T usePartial, Checkerboard *iobj_0)
{
  Checkerboard *board;
  Checkerboard *currentBoard;
  Checkerboard *previousBoard;
  Checkerboard *tmpBoard;
  array<float, 1U> x;
  array<unsigned int, 2U> b_seedIdx;
  array<unsigned int, 2U> seedIdx;
  array<int, 1U> iidx;
  if (scores.size(0) == 0) {
    iobj_0[0].isValid = false;
    iobj_0[0].Energy = rtInfF;
    iobj_0[0].IsDistortionHigh = false;
    board = &iobj_0[0];
    iobj_0[0].BoardIdx.set_size(1, 1);
    iobj_0[0].BoardIdx[0] = 0.0;
    iobj_0[0].BoardIdx.set_size(3, 3);
    for (int i1{0}; i1 < 9; i1++) {
      iobj_0[0].BoardIdx[i1] = 0.0;
    }
    iobj_0[0].BoardCoords.set_size(1, 1, 1);
    iobj_0[0].BoardCoords[0] = 0.0;
    iobj_0[0].BoardCoords.set_size(3, 3, 2);
    for (int i1{0}; i1 < 18; i1++) {
      iobj_0[0].BoardCoords[i1] = 0.0;
    }
    iobj_0[0].Points.set_size(1, 1);
    iobj_0[0].Points[0] = 0.0F;
    iobj_0[0].Points.set_size(0, 2);
  } else {
    int i1;
    int sgn2;
    boolean_T hasExpanded;
    if (points.size(0) < 1) {
      seedIdx.set_size(1, 0);
    } else {
      seedIdx.set_size(1, points.size(0));
      sgn2 = points.size(0) - 1;
      if (static_cast<int>(points.size(0) < 3200)) {
        for (int i{0}; i <= sgn2; i++) {
          seedIdx[i] = static_cast<unsigned int>(i) + 1U;
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int i = 0; i <= sgn2; i++) {
          seedIdx[i] = static_cast<unsigned int>(i) + 1U;
        }
      }
    }
    x.set_size(seedIdx.size(1));
    sgn2 = seedIdx.size(1);
    if (static_cast<int>(seedIdx.size(1) < 3200)) {
      for (int i{0}; i < sgn2; i++) {
        x[i] = scores[static_cast<int>(seedIdx[i]) - 1];
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < sgn2; i++) {
        x[i] = scores[static_cast<int>(seedIdx[i]) - 1];
      }
    }
    ::coder::internal::sort(x, iidx);
    b_seedIdx.set_size(1, iidx.size(0));
    sgn2 = iidx.size(0);
    if (static_cast<int>(iidx.size(0) < 3200)) {
      for (int i{0}; i < sgn2; i++) {
        b_seedIdx[i] = seedIdx[iidx[i] - 1];
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < sgn2; i++) {
        b_seedIdx[i] = seedIdx[iidx[i] - 1];
      }
    }
    seedIdx.set_size(1, b_seedIdx.size(1));
    sgn2 = b_seedIdx.size(1);
    if (static_cast<int>(b_seedIdx.size(1) < 3200)) {
      for (int i{0}; i < sgn2; i++) {
        seedIdx[i] = b_seedIdx[i];
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < sgn2; i++) {
        seedIdx[i] = b_seedIdx[i];
      }
    }
    if (iidx.size(0) > 2000) {
      i1 = static_cast<int>(std::fmin(
          2000.0, std::round(static_cast<double>(seedIdx.size(1)) / 2.0)));
      if (i1 < 1) {
        sgn2 = 0;
      } else {
        sgn2 = i1;
      }
      for (i1 = 0; i1 < sgn2; i1++) {
        seedIdx[i1] = seedIdx[i1];
      }
      seedIdx.set_size(1, sgn2);
    }
    iobj_0[1].isValid = false;
    iobj_0[1].Energy = rtInfF;
    iobj_0[1].IsDistortionHigh = false;
    previousBoard = &iobj_0[1];
    iobj_0[1].BoardIdx.set_size(1, 1);
    iobj_0[1].BoardIdx[0] = 0.0;
    iobj_0[1].BoardIdx.set_size(3, 3);
    for (i1 = 0; i1 < 9; i1++) {
      iobj_0[1].BoardIdx[i1] = 0.0;
    }
    iobj_0[1].BoardCoords.set_size(1, 1, 1);
    iobj_0[1].BoardCoords[0] = 0.0;
    iobj_0[1].BoardCoords.set_size(3, 3, 2);
    for (i1 = 0; i1 < 18; i1++) {
      iobj_0[1].BoardCoords[i1] = 0.0;
    }
    iobj_0[1].Points.set_size(1, 1);
    iobj_0[1].Points[0] = 0.0F;
    iobj_0[1].Points.set_size(0, 2);
    iobj_0[1].IsDistortionHigh = highDistortion;
    iobj_0[2].isValid = false;
    iobj_0[2].Energy = rtInfF;
    iobj_0[2].IsDistortionHigh = false;
    currentBoard = &iobj_0[2];
    iobj_0[2].BoardIdx.set_size(1, 1);
    iobj_0[2].BoardIdx[0] = 0.0;
    iobj_0[2].BoardIdx.set_size(3, 3);
    for (i1 = 0; i1 < 9; i1++) {
      iobj_0[2].BoardIdx[i1] = 0.0;
    }
    iobj_0[2].BoardCoords.set_size(1, 1, 1);
    iobj_0[2].BoardCoords[0] = 0.0;
    iobj_0[2].BoardCoords.set_size(3, 3, 2);
    for (i1 = 0; i1 < 18; i1++) {
      iobj_0[2].BoardCoords[i1] = 0.0;
    }
    iobj_0[2].Points.set_size(1, 1);
    iobj_0[2].Points[0] = 0.0F;
    iobj_0[2].Points.set_size(0, 2);
    iobj_0[2].IsDistortionHigh = highDistortion;
    i1 = seedIdx.size(1);
    for (int b_i{0}; b_i < i1; b_i++) {
      float b_x[2];
      float v2[2];
      float a;
      float ab;
      float adf;
      float b;
      float c;
      float df;
      float tb;
      int i2;
      unsigned int u;
      u = seedIdx[b_i];
      b_x[0] = std::round(points[static_cast<int>(u) - 1]);
      b_x[1] = std::round(points[(static_cast<int>(u) + points.size(0)) - 1]);
      a = Ix2[(static_cast<int>(b_x[1]) +
               Ix2.size(0) * (static_cast<int>(b_x[0]) - 1)) -
              1];
      b = Ixy[(static_cast<int>(b_x[1]) +
               Ixy.size(0) * (static_cast<int>(b_x[0]) - 1)) -
              1];
      c = Iy2[(static_cast<int>(b_x[1]) +
               Iy2.size(0) * (static_cast<int>(b_x[0]) - 1)) -
              1];
      df = a - c;
      adf = std::abs(df);
      tb = b + b;
      ab = std::abs(tb);
      if (adf > ab) {
        b = ab / adf;
        b = adf * std::sqrt(b * b + 1.0F);
      } else if (adf < ab) {
        b = adf / ab;
        b = ab * std::sqrt(b * b + 1.0F);
      } else {
        b = ab * 1.41421354F;
      }
      if (df > 0.0F) {
        adf = df + b;
        sgn2 = 1;
      } else {
        adf = df - b;
        sgn2 = -1;
      }
      if (std::abs(adf) > ab) {
        b = -tb / adf;
        tb = 1.0F / std::sqrt(b * b + 1.0F);
        df = b * tb;
      } else if (ab == 0.0F) {
        df = 1.0F;
        tb = 0.0F;
      } else {
        b = -adf / tb;
        df = 1.0F / std::sqrt(b * b + 1.0F);
        tb = b * df;
      }
      if (a + c < 0.0F) {
        i2 = -1;
      } else {
        i2 = 1;
      }
      if (i2 == sgn2) {
        b = df;
        df = -tb;
        tb = b;
      }
      b_x[0] = -tb * 0.707106769F + df * 0.707106769F;
      adf = -tb * -0.707106769F + df * 0.707106769F;
      b_x[1] = adf;
      v2[0] = df * 0.707106769F + tb * 0.707106769F;
      b = df * -0.707106769F + tb * 0.707106769F;
      v2[1] = b;
      if ((!(std::abs(
                 std::abs(std::abs(rt_atan2f_snf(adf, b_x[0])) - 3.14159274F) -
                 static_cast<float>(theta)) > 0.58904862254808621)) ||
          (!(std::abs(
                 std::abs(std::abs(rt_atan2f_snf(b, v2[0])) - 3.14159274F) -
                 static_cast<float>(theta)) > 0.58904862254808621))) {
        currentBoard->initialize(static_cast<double>(u), points, b_x, v2);
        if (currentBoard->isValid) {
          hasExpanded = true;
          while (hasExpanded) {
            hasExpanded = currentBoard->expandBoardOnce();
          }
        }
        if (currentBoard->Energy < previousBoard->Energy) {
          tmpBoard = previousBoard;
          previousBoard = currentBoard;
          currentBoard = tmpBoard;
        }
      }
    }
    board = previousBoard;
    if (usePartial && previousBoard->isValid) {
      previousBoard->IsDirectionBad[0] = false;
      previousBoard->IsDirectionBad[1] = false;
      previousBoard->IsDirectionBad[2] = false;
      previousBoard->IsDirectionBad[3] = false;
      hasExpanded = true;
      while (hasExpanded) {
        hasExpanded = previousBoard->b_expandBoardOnce();
      }
    }
  }
  return board;
}

//
// Arguments    : Checkerboard *board
//                const ::coder::array<float, 2U> &b_I
// Return Type  : Checkerboard *
//
Checkerboard *orient(Checkerboard *board, const ::coder::array<float, 2U> &b_I)
{
  Checkerboard *b_board;
  array<double, 3U> c_y;
  array<double, 3U> r;
  array<double, 3U> r1;
  array<double, 2U> c_board;
  array<double, 2U> newBoardCoords1;
  array<double, 2U> newBoardCoords2;
  array<boolean_T, 3U> b_x;
  array<boolean_T, 2U> y;
  float x;
  int numRot_size[2];
  b_board = board;
  x = b_board->Energy;
  if (!std::isinf(x)) {
    double numRot_data;
    int i;
    int i2;
    int ix;
    int iy;
    int loop_ub;
    boolean_T b_y;
    boolean_T exitg1;
    if (b_board->BoardCoords.size(0) < b_board->BoardCoords.size(1)) {
      c_board.set_size(b_board->BoardIdx.size(0), b_board->BoardIdx.size(1));
      i2 = b_board->BoardIdx.size(0) * b_board->BoardIdx.size(1) - 1;
      for (i = 0; i <= i2; i++) {
        c_board[i] = b_board->BoardIdx[i];
      }
      rot90(c_board, b_board->BoardIdx);
      c_board.set_size(b_board->BoardCoords.size(0),
                       b_board->BoardCoords.size(1));
      i2 = b_board->BoardCoords.size(1);
      for (i = 0; i < i2; i++) {
        loop_ub = b_board->BoardCoords.size(0);
        for (ix = 0; ix < loop_ub; ix++) {
          c_board[ix + c_board.size(0) * i] =
              b_board->BoardCoords[ix + b_board->BoardCoords.size(0) * i];
        }
      }
      rot90(c_board, newBoardCoords1);
      c_board.set_size(b_board->BoardCoords.size(0),
                       b_board->BoardCoords.size(1));
      i2 = b_board->BoardCoords.size(1);
      for (i = 0; i < i2; i++) {
        loop_ub = b_board->BoardCoords.size(0);
        for (ix = 0; ix < loop_ub; ix++) {
          c_board[ix + c_board.size(0) * i] =
              b_board->BoardCoords[(ix + b_board->BoardCoords.size(0) * i) +
                                   b_board->BoardCoords.size(0) *
                                       b_board->BoardCoords.size(1)];
        }
      }
      rot90(c_board, newBoardCoords2);
      c_y.set_size(newBoardCoords1.size(0), newBoardCoords1.size(1), 2);
      i = newBoardCoords1.size(0) * newBoardCoords1.size(1);
      if (static_cast<int>(i < 3200)) {
        for (int j{0}; j < i; j++) {
          c_y[j] = newBoardCoords1[j];
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int j = 0; j < i; j++) {
          c_y[j] = newBoardCoords1[j];
        }
      }
      iy = i + -1;
      ix = newBoardCoords2.size(0) * newBoardCoords2.size(1);
      if (static_cast<int>(ix < 3200)) {
        for (int j{0}; j < ix; j++) {
          c_y[i + j] = newBoardCoords2[j];
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int j = 0; j < ix; j++) {
          c_y[(iy + j) + 1] = newBoardCoords2[j];
        }
      }
      b_board->BoardCoords.set_size(c_y.size(0), c_y.size(1), 2);
      i2 = c_y.size(0) * c_y.size(1) * 2;
      for (i = 0; i < i2; i++) {
        b_board->BoardCoords[i] = c_y[i];
      }
    }
    newBoardCoords1.set_size(b_board->BoardIdx.size(0),
                             b_board->BoardIdx.size(1));
    i2 = b_board->BoardIdx.size(0) * b_board->BoardIdx.size(1);
    for (i = 0; i < i2; i++) {
      newBoardCoords1[i] = b_board->BoardIdx[i];
    }
    y.set_size(1, newBoardCoords1.size(1));
    i2 = newBoardCoords1.size(1);
    if (static_cast<int>(newBoardCoords1.size(1) < 3200)) {
      for (int j{0}; j < i2; j++) {
        y[j] = true;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int j = 0; j < i2; j++) {
        y[j] = true;
      }
    }
    iy = newBoardCoords1.size(1);
    i2 = 0;
    for (i = 0; i < iy; i++) {
      loop_ub = i2 + newBoardCoords1.size(0);
      ix = i2;
      i2 += newBoardCoords1.size(0);
      exitg1 = false;
      while ((!exitg1) && (ix + 1 <= loop_ub)) {
        if (newBoardCoords1[ix] == 0.0) {
          y[i] = false;
          exitg1 = true;
        } else {
          ix++;
        }
      }
    }
    b_y = (y.size(1) != 0);
    if (b_y) {
      iy = 0;
      exitg1 = false;
      while ((!exitg1) && (iy <= y.size(1) - 1)) {
        if (!y[iy]) {
          b_y = false;
          exitg1 = true;
        } else {
          iy++;
        }
      }
    }
    if (b_y) {
      if (!isUpperLeftBlack(b_board, b_I)) {
        c_board.set_size(b_board->BoardIdx.size(0), b_board->BoardIdx.size(1));
        i2 = b_board->BoardIdx.size(0) * b_board->BoardIdx.size(1) - 1;
        for (i = 0; i <= i2; i++) {
          c_board[i] = b_board->BoardIdx[i];
        }
        b_rot90(c_board, b_board->BoardIdx);
        c_board.set_size(b_board->BoardCoords.size(0),
                         b_board->BoardCoords.size(1));
        i2 = b_board->BoardCoords.size(1);
        for (i = 0; i < i2; i++) {
          loop_ub = b_board->BoardCoords.size(0);
          for (ix = 0; ix < loop_ub; ix++) {
            c_board[ix + c_board.size(0) * i] =
                b_board->BoardCoords[ix + b_board->BoardCoords.size(0) * i];
          }
        }
        b_rot90(c_board, newBoardCoords1);
        c_board.set_size(b_board->BoardCoords.size(0),
                         b_board->BoardCoords.size(1));
        i2 = b_board->BoardCoords.size(1);
        for (i = 0; i < i2; i++) {
          loop_ub = b_board->BoardCoords.size(0);
          for (ix = 0; ix < loop_ub; ix++) {
            c_board[ix + c_board.size(0) * i] =
                b_board->BoardCoords[(ix + b_board->BoardCoords.size(0) * i) +
                                     b_board->BoardCoords.size(0) *
                                         b_board->BoardCoords.size(1)];
          }
        }
        b_rot90(c_board, newBoardCoords2);
        c_y.set_size(newBoardCoords1.size(0), newBoardCoords1.size(1), 2);
        i = newBoardCoords1.size(0) * newBoardCoords1.size(1);
        if (static_cast<int>(i < 3200)) {
          for (int j{0}; j < i; j++) {
            c_y[j] = newBoardCoords1[j];
          }
        } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

          for (int j = 0; j < i; j++) {
            c_y[j] = newBoardCoords1[j];
          }
        }
        iy = i + -1;
        ix = newBoardCoords2.size(0) * newBoardCoords2.size(1);
        if (static_cast<int>(ix < 3200)) {
          for (int j{0}; j < ix; j++) {
            c_y[i + j] = newBoardCoords2[j];
          }
        } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

          for (int j = 0; j < ix; j++) {
            c_y[(iy + j) + 1] = newBoardCoords2[j];
          }
        }
        b_board->BoardCoords.set_size(c_y.size(0), c_y.size(1), 2);
        i2 = c_y.size(0) * c_y.size(1) * 2;
        for (i = 0; i < i2; i++) {
          b_board->BoardCoords[i] = c_y[i];
        }
      }
    } else {
      double cornerIdx[4];
      signed char ii_data;
      cornerIdx[0] = b_board->BoardIdx[0];
      cornerIdx[1] = b_board->BoardIdx[b_board->BoardIdx.size(0) - 1];
      cornerIdx[2] = b_board->BoardIdx[(b_board->BoardIdx.size(0) +
                                        b_board->BoardIdx.size(0) *
                                            (b_board->BoardIdx.size(1) - 1)) -
                                       1];
      cornerIdx[3] = b_board->BoardIdx[b_board->BoardIdx.size(0) *
                                       (b_board->BoardIdx.size(1) - 1)];
      iy = 0;
      i2 = 1;
      loop_ub = 0;
      exitg1 = false;
      while ((!exitg1) && (loop_ub < 4)) {
        if (cornerIdx[loop_ub] != 0.0) {
          iy = 1;
          ii_data = static_cast<signed char>(loop_ub + 1);
          exitg1 = true;
        } else {
          loop_ub++;
        }
      }
      if (iy == 0) {
        i2 = 0;
      }
      numRot_size[0] = 1;
      numRot_size[1] = i2;
      for (i = 0; i < i2; i++) {
        numRot_data = static_cast<double>(ii_data) - 1.0;
      }
      if (i2 == 0) {
        numRot_size[0] = 1;
        numRot_size[1] = 1;
        numRot_data = 0.0;
      }
      y.set_size(1, 1);
      y[0] = (numRot_data == 2.0);
      if (!y[0]) {
        if (isequal((const double *)&numRot_data, numRot_size, 1.0)) {
          iy = 1;
        } else if (isequal((const double *)&numRot_data, numRot_size, 3.0)) {
          iy = 2;
        } else {
          iy = 0;
        }
        if (iy != 0) {
          newBoardCoords1.set_size(b_board->BoardIdx.size(0),
                                   b_board->BoardIdx.size(1));
          i2 = b_board->BoardIdx.size(0) * b_board->BoardIdx.size(1);
          for (i = 0; i < i2; i++) {
            newBoardCoords1[i] = b_board->BoardIdx[i];
          }
          flip(newBoardCoords1, static_cast<double>(iy));
          b_board->BoardIdx.set_size(newBoardCoords1.size(0),
                                     newBoardCoords1.size(1));
          i2 = newBoardCoords1.size(0) * newBoardCoords1.size(1);
          for (i = 0; i < i2; i++) {
            b_board->BoardIdx[i] = newBoardCoords1[i];
          }
          newBoardCoords1.set_size(b_board->BoardCoords.size(0),
                                   b_board->BoardCoords.size(1));
          i2 = b_board->BoardCoords.size(1);
          for (i = 0; i < i2; i++) {
            loop_ub = b_board->BoardCoords.size(0);
            for (ix = 0; ix < loop_ub; ix++) {
              newBoardCoords1[ix + newBoardCoords1.size(0) * i] =
                  b_board->BoardCoords[ix + b_board->BoardCoords.size(0) * i];
            }
          }
          flip(newBoardCoords1, static_cast<double>(iy));
          newBoardCoords2.set_size(b_board->BoardCoords.size(0),
                                   b_board->BoardCoords.size(1));
          i2 = b_board->BoardCoords.size(1);
          for (i = 0; i < i2; i++) {
            loop_ub = b_board->BoardCoords.size(0);
            for (ix = 0; ix < loop_ub; ix++) {
              newBoardCoords2[ix + newBoardCoords2.size(0) * i] =
                  b_board->BoardCoords[(ix + b_board->BoardCoords.size(0) * i) +
                                       b_board->BoardCoords.size(0) *
                                           b_board->BoardCoords.size(1)];
            }
          }
          flip(newBoardCoords2, static_cast<double>(iy));
          c_y.set_size(newBoardCoords1.size(0), newBoardCoords1.size(1), 2);
          i = newBoardCoords1.size(0) * newBoardCoords1.size(1);
          if (static_cast<int>(i < 3200)) {
            for (int j{0}; j < i; j++) {
              c_y[j] = newBoardCoords1[j];
            }
          } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

            for (int j = 0; j < i; j++) {
              c_y[j] = newBoardCoords1[j];
            }
          }
          iy = i + -1;
          ix = newBoardCoords2.size(0) * newBoardCoords2.size(1);
          if (static_cast<int>(ix < 3200)) {
            for (int j{0}; j < ix; j++) {
              c_y[i + j] = newBoardCoords2[j];
            }
          } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

            for (int j = 0; j < ix; j++) {
              c_y[(iy + j) + 1] = newBoardCoords2[j];
            }
          }
          b_board->BoardCoords.set_size(c_y.size(0), c_y.size(1), 2);
          i2 = c_y.size(0) * c_y.size(1) * 2;
          for (i = 0; i < i2; i++) {
            b_board->BoardCoords[i] = c_y[i];
          }
        }
      } else {
        c_board.set_size(b_board->BoardIdx.size(0), b_board->BoardIdx.size(1));
        i2 = b_board->BoardIdx.size(0) * b_board->BoardIdx.size(1) - 1;
        for (i = 0; i <= i2; i++) {
          c_board[i] = b_board->BoardIdx[i];
        }
        b_rot90(c_board, b_board->BoardIdx);
        c_board.set_size(b_board->BoardCoords.size(0),
                         b_board->BoardCoords.size(1));
        i2 = b_board->BoardCoords.size(1);
        for (i = 0; i < i2; i++) {
          loop_ub = b_board->BoardCoords.size(0);
          for (ix = 0; ix < loop_ub; ix++) {
            c_board[ix + c_board.size(0) * i] =
                b_board->BoardCoords[ix + b_board->BoardCoords.size(0) * i];
          }
        }
        b_rot90(c_board, newBoardCoords1);
        c_board.set_size(b_board->BoardCoords.size(0),
                         b_board->BoardCoords.size(1));
        i2 = b_board->BoardCoords.size(1);
        for (i = 0; i < i2; i++) {
          loop_ub = b_board->BoardCoords.size(0);
          for (ix = 0; ix < loop_ub; ix++) {
            c_board[ix + c_board.size(0) * i] =
                b_board->BoardCoords[(ix + b_board->BoardCoords.size(0) * i) +
                                     b_board->BoardCoords.size(0) *
                                         b_board->BoardCoords.size(1)];
          }
        }
        b_rot90(c_board, newBoardCoords2);
        c_y.set_size(newBoardCoords1.size(0), newBoardCoords1.size(1), 2);
        i = newBoardCoords1.size(0) * newBoardCoords1.size(1);
        if (static_cast<int>(i < 3200)) {
          for (int j{0}; j < i; j++) {
            c_y[j] = newBoardCoords1[j];
          }
        } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

          for (int j = 0; j < i; j++) {
            c_y[j] = newBoardCoords1[j];
          }
        }
        iy = i + -1;
        ix = newBoardCoords2.size(0) * newBoardCoords2.size(1);
        if (static_cast<int>(ix < 3200)) {
          for (int j{0}; j < ix; j++) {
            c_y[i + j] = newBoardCoords2[j];
          }
        } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

          for (int j = 0; j < ix; j++) {
            c_y[(iy + j) + 1] = newBoardCoords2[j];
          }
        }
        b_board->BoardCoords.set_size(c_y.size(0), c_y.size(1), 2);
        i2 = c_y.size(0) * c_y.size(1) * 2;
        for (i = 0; i < i2; i++) {
          b_board->BoardCoords[i] = c_y[i];
        }
      }
    }
    iy = b_board->BoardCoords.size(0);
    loop_ub = b_board->BoardCoords.size(1);
    if (iy == 0) {
      iy = 0;
    } else {
      iy = static_cast<int>(std::fmod(static_cast<double>(iy), 2.0));
    }
    if (loop_ub == 0) {
      i2 = 0;
    } else {
      i2 = static_cast<int>(std::fmod(static_cast<double>(loop_ub), 2.0));
    }
    if ((iy == 0.0) == (i2 == 0.0)) {
      r.set_size(1, 1, b_board->BoardCoords.size(2));
      i2 = b_board->BoardCoords.size(2);
      for (i = 0; i < i2; i++) {
        r[i] = b_board->BoardCoords[b_board->BoardCoords.size(0) *
                                    b_board->BoardCoords.size(1) * i];
      }
      loop_ub = b_board->BoardCoords.size(0);
      iy = b_board->BoardCoords.size(1);
      r1.set_size(1, 1, b_board->BoardCoords.size(2));
      i2 = b_board->BoardCoords.size(2);
      for (i = 0; i < i2; i++) {
        r1[i] = b_board->BoardCoords[((loop_ub + b_board->BoardCoords.size(0) *
                                                     (iy - 1)) +
                                      b_board->BoardCoords.size(0) *
                                          b_board->BoardCoords.size(1) * i) -
                                     1];
      }
      b_x.set_size(1, 1, r.size(2));
      i2 = r.size(2);
      if (static_cast<int>(r.size(2) < 3200)) {
        for (int j{0}; j < i2; j++) {
          b_x[j] = (r[j] > r1[j]);
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int j = 0; j < i2; j++) {
          b_x[j] = (r[j] > r1[j]);
        }
      }
      b_y = false;
      ix = 1;
      exitg1 = false;
      while ((!exitg1) && (ix <= b_x.size(2))) {
        if (b_x[ix - 1]) {
          b_y = true;
          exitg1 = true;
        } else {
          ix++;
        }
      }
      if (b_y) {
        numRot_data =
            b_board->BoardCoords[(b_board->BoardCoords.size(0) +
                                  b_board->BoardCoords.size(0) *
                                      (b_board->BoardCoords.size(1) - 1)) -
                                 1];
        if (numRot_data != 0.0) {
          c_board.set_size(b_board->BoardIdx.size(0),
                           b_board->BoardIdx.size(1));
          i2 = b_board->BoardIdx.size(0) * b_board->BoardIdx.size(1) - 1;
          for (i = 0; i <= i2; i++) {
            c_board[i] = b_board->BoardIdx[i];
          }
          b_rot90(c_board, b_board->BoardIdx);
          c_board.set_size(b_board->BoardCoords.size(0),
                           b_board->BoardCoords.size(1));
          i2 = b_board->BoardCoords.size(1);
          for (i = 0; i < i2; i++) {
            loop_ub = b_board->BoardCoords.size(0);
            for (ix = 0; ix < loop_ub; ix++) {
              c_board[ix + c_board.size(0) * i] =
                  b_board->BoardCoords[ix + b_board->BoardCoords.size(0) * i];
            }
          }
          b_rot90(c_board, newBoardCoords1);
          c_board.set_size(b_board->BoardCoords.size(0),
                           b_board->BoardCoords.size(1));
          i2 = b_board->BoardCoords.size(1);
          for (i = 0; i < i2; i++) {
            loop_ub = b_board->BoardCoords.size(0);
            for (ix = 0; ix < loop_ub; ix++) {
              c_board[ix + c_board.size(0) * i] =
                  b_board->BoardCoords[(ix + b_board->BoardCoords.size(0) * i) +
                                       b_board->BoardCoords.size(0) *
                                           b_board->BoardCoords.size(1)];
            }
          }
          b_rot90(c_board, newBoardCoords2);
          c_y.set_size(newBoardCoords1.size(0), newBoardCoords1.size(1), 2);
          i = newBoardCoords1.size(0) * newBoardCoords1.size(1);
          if (static_cast<int>(i < 3200)) {
            for (int j{0}; j < i; j++) {
              c_y[j] = newBoardCoords1[j];
            }
          } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

            for (int j = 0; j < i; j++) {
              c_y[j] = newBoardCoords1[j];
            }
          }
          iy = i + -1;
          ix = newBoardCoords2.size(0) * newBoardCoords2.size(1);
          if (static_cast<int>(ix < 3200)) {
            for (int j{0}; j < ix; j++) {
              c_y[i + j] = newBoardCoords2[j];
            }
          } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

            for (int j = 0; j < ix; j++) {
              c_y[(iy + j) + 1] = newBoardCoords2[j];
            }
          }
          b_board->BoardCoords.set_size(c_y.size(0), c_y.size(1), 2);
          i2 = c_y.size(0) * c_y.size(1) * 2;
          for (i = 0; i < i2; i++) {
            b_board->BoardCoords[i] = c_y[i];
          }
        }
      }
    }
  }
  return b_board;
}

//
// Arguments    : const Checkerboard *b_this
//                boolean_T usePartial
//                ::coder::array<double, 2U> &points
//                double boardSize[2]
// Return Type  : void
//
void toPoints(const Checkerboard *b_this, boolean_T usePartial,
              ::coder::array<double, 2U> &points, double boardSize[2])
{
  array<double, 2U> b_x;
  array<boolean_T, 1U> x;
  int b_loop_ub;
  int i1;
  int loop_ub;
  boolean_T exitg1;
  boolean_T y;
  x.set_size(b_this->BoardIdx.size(0) * b_this->BoardIdx.size(1));
  loop_ub = b_this->BoardIdx.size(0) * b_this->BoardIdx.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      x[i] = (b_this->BoardIdx[i] == 0.0);
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      x[i] = (b_this->BoardIdx[i] == 0.0);
    }
  }
  y = false;
  loop_ub = 1;
  exitg1 = false;
  while ((!exitg1) && (loop_ub <= x.size(0))) {
    if (x[loop_ub - 1]) {
      y = true;
      exitg1 = true;
    } else {
      loop_ub++;
    }
  }
  if (y && (!usePartial)) {
    points.set_size(0, 0);
    boardSize[0] = 0.0;
    boardSize[1] = 0.0;
  } else {
    double numPoints;
    numPoints = static_cast<double>(b_this->BoardCoords.size(0)) *
                static_cast<double>(b_this->BoardCoords.size(1));
    points.set_size(static_cast<int>(numPoints), 2);
    loop_ub = static_cast<int>(numPoints) << 1;
    if (static_cast<int>(loop_ub < 3200)) {
      for (int i{0}; i < loop_ub; i++) {
        points[i] = 0.0;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < loop_ub; i++) {
        points[i] = 0.0;
      }
    }
    b_x.set_size(b_this->BoardCoords.size(1), b_this->BoardCoords.size(0));
    loop_ub = b_this->BoardCoords.size(0);
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(i1, b_loop_ub)

    for (int i = 0; i < loop_ub; i++) {
      b_loop_ub = b_this->BoardCoords.size(1);
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        b_x[i1 + b_x.size(0) * i] =
            b_this->BoardCoords[i + b_this->BoardCoords.size(0) * i1];
      }
    }
    loop_ub = b_x.size(0) * b_x.size(1);
    if (static_cast<int>(loop_ub < 3200)) {
      for (int i{0}; i < loop_ub; i++) {
        points[i] = b_x[i];
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < loop_ub; i++) {
        points[i] = b_x[i];
      }
    }
    b_x.set_size(b_this->BoardCoords.size(1), b_this->BoardCoords.size(0));
    loop_ub = b_this->BoardCoords.size(0);
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(i1, b_loop_ub)

    for (int i = 0; i < loop_ub; i++) {
      b_loop_ub = b_this->BoardCoords.size(1);
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        b_x[i1 + b_x.size(0) * i] =
            b_this->BoardCoords[(i + b_this->BoardCoords.size(0) * i1) +
                                b_this->BoardCoords.size(0) *
                                    b_this->BoardCoords.size(1)];
      }
    }
    loop_ub = b_x.size(0) * b_x.size(1);
    if (static_cast<int>(loop_ub < 3200)) {
      for (int i{0}; i < loop_ub; i++) {
        points[i + points.size(0)] = b_x[i];
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < loop_ub; i++) {
        points[i + points.size(0)] = b_x[i];
      }
    }
    boardSize[0] = static_cast<unsigned int>(b_this->BoardCoords.size(1)) + 1U;
    boardSize[1] = static_cast<unsigned int>(b_this->BoardCoords.size(0)) + 1U;
  }
}

} // namespace checkerboard
} // namespace calibration
} // namespace internal
} // namespace vision
} // namespace coder

//
// File trailer for detectCheckerboard.cpp
//
// [EOF]
//
