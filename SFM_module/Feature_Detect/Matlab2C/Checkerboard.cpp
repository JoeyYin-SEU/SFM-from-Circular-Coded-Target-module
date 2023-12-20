//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: Checkerboard.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "Checkerboard.h"
#include "binaryMinOrMax.h"
#include "bsxfun.h"
#include "colon.h"
#include "combineVectorElements.h"
#include "div.h"
#include "eml_setop.h"
#include "find.h"
#include "get_chessborad_pixel_rtwutil.h"
#include "ismember.h"
#include "minOrMax.h"
#include "norm.h"
#include "polyfit.h"
#include "polyval.h"
#include "rt_nonfinite.h"
#include "squeeze.h"
#include "unsafeSxfun.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>

// Function Declarations
static void b_binary_expand_op(coder::vision::internal::calibration::
  checkerboard::Checkerboard *in1, const coder::array<float, 2U> &in2, const
  coder::array<float, 1U> &in3, const coder::array<float, 2U> &in4, const coder::
  array<float, 2U> &in5);
static void b_minus(coder::array<float, 2U> &in1, const coder::array<float, 2U>
                    &in2);
static void binary_expand_op(coder::vision::internal::calibration::checkerboard::
  Checkerboard *in1, const coder::array<float, 2U> &in2, const coder::array<
  float, 1U> &in3, const coder::array<float, 2U> &in4, const coder::array<float,
  2U> &in5);
static void binary_expand_op(coder::array<float, 1U> &in1, const coder::array<
  float, 1U> &in2, const coder::array<float, 1U> &in3);
static void binary_expand_op(coder::array<double, 2U> &in1, const coder::vision::
  internal::calibration::checkerboard::Checkerboard *in2, const coder::array<int,
  1U> &in3, const coder::array<int, 1U> &in4, const coder::array<int, 1U> &in5);
static void binary_expand_op(coder::array<double, 2U> &in1, const coder::array<
  double, 2U> &in2, const coder::array<double, 2U> &in3);
static void binary_expand_op(coder::array<double, 2U> &in1, const coder::array<
  double, 3U> &in2, const coder::array<double, 3U> &in3, const coder::array<
  double, 3U> &in4);
static void binary_expand_op(coder::array<double, 2U> &in1, const coder::vision::
  internal::calibration::checkerboard::Checkerboard *in2, const coder::array<
  double, 2U> &in3, const coder::array<double, 2U> &in4);
static void binary_expand_op(coder::array<double, 2U> &in1, const coder::vision::
  internal::calibration::checkerboard::Checkerboard *in2, const coder::array<int,
  1U> &in3, const coder::array<double, 2U> &in4, const coder::array<int, 1U>
  &in5, const coder::array<int, 1U> &in6);
static void c_binary_expand_op(coder::vision::internal::calibration::
  checkerboard::Checkerboard *in1, const coder::array<float, 2U> &in2, const
  coder::array<float, 1U> &in3, const coder::array<float, 2U> &in4, const coder::
  array<float, 2U> &in5);
static void c_binary_expand_op(coder::array<float, 2U> &in1, const coder::array<
  float, 2U> &in2, const coder::array<float, 2U> &in3);
static void d_binary_expand_op(coder::vision::internal::calibration::
  checkerboard::Checkerboard *in1, const coder::array<float, 2U> &in2, const
  coder::array<float, 1U> &in3, const coder::array<float, 2U> &in4, const coder::
  array<float, 2U> &in5);
static void e_binary_expand_op(coder::array<boolean_T, 1U> &in1, const coder::
  array<double, 2U> &in2, const coder::array<double, 2U> &in3);
static void minus(coder::array<float, 2U> &in1, const coder::array<float, 2U>
                  &in2);

// Function Definitions
//
// Arguments    : const ::coder::array<boolean_T, 2U> &arr
//                ::coder::array<double, 2U> &matchedIdx
// Return Type  : void
//
namespace coder
{
  namespace vision
  {
    namespace internal
    {
      namespace calibration
      {
        namespace checkerboard
        {
          void Checkerboard::arrayFind(const ::coder::array<boolean_T, 2U> &arr,
            ::coder::array<double, 2U> &matchedIdx)
          {
            array<int, 2U> b_ii;
            array<signed char, 2U> matchArr;
            array<boolean_T, 2U> x;
            int idx;
            int ii;
            int nx;
            boolean_T exitg1;
            matchArr.set_size(1, arr.size(1) - 2);
            ii = arr.size(1);
            for (idx = 0; idx <= ii - 3; idx++) {
              boolean_T b_x[3];
              boolean_T y;
              for (nx = 0; nx < 3; nx++) {
                b_x[nx] = arr[idx + nx];
              }

              y = true;
              nx = 0;
              exitg1 = false;
              while ((!exitg1) && (nx <= 2)) {
                if (!b_x[nx]) {
                  y = false;
                  exitg1 = true;
                } else {
                  nx++;
                }
              }

              matchArr[idx] = static_cast<signed char>(y);
            }

            x.set_size(1, matchArr.size(1));
            nx = matchArr.size(1);
            if (static_cast<int>(matchArr.size(1) < 3200)) {
              for (int i{0}; i < nx; i++) {
                x[i] = (matchArr[i] == 1);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < nx; i++) {
                x[i] = (matchArr[i] == 1);
              }
            }

            nx = x.size(1);
            idx = 0;
            b_ii.set_size(1, x.size(1));
            ii = 0;
            exitg1 = false;
            while ((!exitg1) && (ii <= nx - 1)) {
              if (x[ii]) {
                idx++;
                b_ii[idx - 1] = ii + 1;
                if (idx >= nx) {
                  exitg1 = true;
                } else {
                  ii++;
                }
              } else {
                ii++;
              }
            }

            if (x.size(1) == 1) {
              if (idx == 0) {
                b_ii.set_size(1, 0);
              }
            } else {
              if (idx < 1) {
                idx = 0;
              }

              b_ii.set_size(b_ii.size(0), idx);
            }

            matchedIdx.set_size(1, b_ii.size(1));
            nx = b_ii.size(1);
            if (static_cast<int>(b_ii.size(1) < 3200)) {
              for (int i{0}; i < nx; i++) {
                matchedIdx[i] = b_ii[i];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < nx; i++) {
                matchedIdx[i] = b_ii[i];
              }
            }
          }

          //
          // Arguments    : const ::coder::array<double, 2U> &predictedPoints
          //                ::coder::array<double, 2U> &indices
          // Return Type  : void
          //
          void Checkerboard::b_findClosestIndices(const ::coder::array<double,
            2U> &predictedPoints, ::coder::array<double, 2U> &indices) const
          {
            array<double, 2U> b_indices;
            array<double, 2U> p;
            array<double, 2U> remIdx;
            array<double, 1U> b_this;
            array<double, 1U> validBoardIdx;
            array<float, 2U> c_this;
            array<float, 2U> diffs;
            array<float, 1U> dists;
            array<int, 2U> r1;
            array<int, 1U> validPredictions;
            array<boolean_T, 2U> distIdx;
            array<boolean_T, 1U> r;
            float minDist;
            int b_loop_ub;
            int c_loop_ub;
            int d_loop_ub;
            int i1;
            int iindx;
            int loop_ub;
            int n;
            int nx;
            indices.set_size(1, predictedPoints.size(0));
            loop_ub = predictedPoints.size(0);
            if (static_cast<int>(predictedPoints.size(0) < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                indices[i] = 0.0;
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                indices[i] = 0.0;
              }
            }

            nx = Points.size(0);
            if (nx < 1) {
              p.set_size(1, 0);
            } else {
              p.set_size(1, nx);
              loop_ub = nx - 1;
              if (static_cast<int>(nx < 3200)) {
                for (int i{0}; i <= loop_ub; i++) {
                  p[i] = static_cast<double>(i) + 1.0;
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i <= loop_ub; i++) {
                  p[i] = static_cast<double>(i) + 1.0;
                }
              }
            }

            nx = BoardIdx.size(0) * BoardIdx.size(1);
            b_this = BoardIdx.reshape(nx);
            do_vectors(p, b_this, remIdx, validPredictions, &nx);
            if (remIdx.size(1) != 0) {
              int b_i;
              r.set_size(predictedPoints.size(0));
              loop_ub = predictedPoints.size(0);
              if (static_cast<int>(predictedPoints.size(0) < 3200)) {
                for (int i{0}; i < loop_ub; i++) {
                  r[i] = std::isnan(predictedPoints[i]);
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < loop_ub; i++) {
                  r[i] = std::isnan(predictedPoints[i]);
                }
              }

              loop_ub = r.size(0);
              if (static_cast<int>(r.size(0) < 3200)) {
                for (int i{0}; i < loop_ub; i++) {
                  r[i] = !r[i];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < loop_ub; i++) {
                  r[i] = !r[i];
                }
              }

              b_eml_find(r, validPredictions);
              b_i = validPredictions.size(0);
              if (validPredictions.size(0) - 1 >= 0) {
                b_loop_ub = predictedPoints.size(1);
                c_loop_ub = Points.size(1);
                n = BoardIdx.size(0) * BoardIdx.size(1);
                i1 = BoardIdx.size(0) * BoardIdx.size(1);
                d_loop_ub = Points.size(1);
              }

              for (int c_i{0}; c_i < b_i; c_i++) {
                int d_i;
                p.set_size(1, predictedPoints.size(1));
                for (nx = 0; nx < b_loop_ub; nx++) {
                  p[nx] = predictedPoints[(validPredictions[c_i] +
                    predictedPoints.size(0) * nx) - 1];
                }

                c_this.set_size(remIdx.size(1), Points.size(1));
                for (nx = 0; nx < c_loop_ub; nx++) {
                  loop_ub = remIdx.size(1);
                  for (d_i = 0; d_i < loop_ub; d_i++) {
                    c_this[d_i + c_this.size(0) * nx] = Points[(static_cast<int>
                      (remIdx[d_i]) + Points.size(0) * nx) - 1];
                  }
                }

                bsxfun(c_this, p, diffs);
                dists.set_size(diffs.size(0));
                nx = diffs.size(0);
                for (loop_ub = 0; loop_ub < nx; loop_ub++) {
                  dists[loop_ub] = rt_hypotf_snf(diffs[loop_ub], diffs[loop_ub +
                    diffs.size(0)]);
                }

                loop_ub = indices.size(1) - 1;
                nx = 0;
                for (d_i = 0; d_i <= loop_ub; d_i++) {
                  if (indices[d_i] > 0.0) {
                    nx++;
                  }
                }

                r1.set_size(1, nx);
                nx = 0;
                for (d_i = 0; d_i <= loop_ub; d_i++) {
                  if (indices[d_i] > 0.0) {
                    r1[nx] = d_i + 1;
                    nx++;
                  }
                }

                b_indices.set_size(1, r1.size(1));
                loop_ub = r1.size(1);
                for (nx = 0; nx < loop_ub; nx++) {
                  b_indices[nx] = indices[r1[nx] - 1];
                }

                isMember(remIdx, b_indices, distIdx);
                loop_ub = distIdx.size(1);
                for (d_i = 0; d_i < loop_ub; d_i++) {
                  if (distIdx[d_i]) {
                    dists[d_i] = rtInfF;
                  }
                }

                ::coder::internal::minimum(dists, &minDist, &iindx);
                nx = 0;
                for (loop_ub = 0; loop_ub < i1; loop_ub++) {
                  if (BoardIdx[loop_ub] != 0.0) {
                    nx++;
                  }
                }

                validBoardIdx.set_size(nx);
                d_i = -1;
                for (loop_ub = 0; loop_ub < n; loop_ub++) {
                  if (BoardIdx[loop_ub] != 0.0) {
                    d_i++;
                    validBoardIdx[d_i] = BoardIdx[loop_ub];
                  }
                }

                c_this.set_size(validBoardIdx.size(0), Points.size(1));
                for (nx = 0; nx < d_loop_ub; nx++) {
                  loop_ub = validBoardIdx.size(0);
                  for (d_i = 0; d_i < loop_ub; d_i++) {
                    c_this[d_i + c_this.size(0) * nx] = Points[(static_cast<int>
                      (validBoardIdx[d_i]) + Points.size(0) * nx) - 1];
                  }
                }

                bsxfun(c_this, p, diffs);
                dists.set_size(diffs.size(0));
                nx = diffs.size(0);
                for (loop_ub = 0; loop_ub < nx; loop_ub++) {
                  dists[loop_ub] = rt_hypotf_snf(diffs[loop_ub], diffs[loop_ub +
                    diffs.size(0)]);
                }

                if (minDist < ::coder::internal::minimum(dists) / 2.0F) {
                  indices[validPredictions[c_i] - 1] = remIdx[iindx - 1];
                }
              }

              n = 0;
              b_i = indices.size(1);
              for (loop_ub = 0; loop_ub < b_i; loop_ub++) {
                if (indices[loop_ub] != 0.0) {
                  n++;
                }
              }

              if (n < 4) {
                loop_ub = indices.size(1);
                if (static_cast<int>(indices.size(1) < 3200)) {
                  for (int i{0}; i < loop_ub; i++) {
                    if (indices[i] > 0.0) {
                      indices[i] = 0.0;
                    }
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (int i = 0; i < loop_ub; i++) {
                    if (indices[i] > 0.0) {
                      indices[i] = 0.0;
                    }
                  }
                }
              }
            }
          }

          //
          // Arguments    : double coordsToUse[2]
          // Return Type  : void
          //
          void Checkerboard::b_findIndependentVar(double coordsToUse[2]) const
          {
            array<double, 1U> b_x;
            array<double, 1U> x;
            array<int, 1U> r2;
            array<int, 1U> r3;
            array<int, 1U> r4;
            array<int, 1U> r5;
            array<boolean_T, 1U> r;
            array<boolean_T, 1U> r1;
            int end;
            int loop_ub;
            r.set_size(BoardIdx.size(0));
            loop_ub = BoardIdx.size(0);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                r[i] = (BoardIdx[i] > 0.0);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                r[i] = (BoardIdx[i] > 0.0);
              }
            }

            r1.set_size(BoardIdx.size(0));
            loop_ub = BoardIdx.size(0);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                r1[i] = (BoardIdx[i + BoardIdx.size(0)] > 0.0);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                r1[i] = (BoardIdx[i + BoardIdx.size(0)] > 0.0);
              }
            }

            end = r.size(0) - 1;
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (r[b_i] && r1[b_i]) {
                loop_ub++;
              }
            }

            r2.set_size(loop_ub);
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (r[b_i] && r1[b_i]) {
                r2[loop_ub] = b_i + 1;
                loop_ub++;
              }
            }

            end = r.size(0) - 1;
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (r[b_i] && r1[b_i]) {
                loop_ub++;
              }
            }

            r3.set_size(loop_ub);
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (r[b_i] && r1[b_i]) {
                r3[loop_ub] = b_i + 1;
                loop_ub++;
              }
            }

            x.set_size(r2.size(0));
            loop_ub = r2.size(0);
            if (static_cast<int>(r2.size(0) < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                x[i] = BoardCoords[(r2[i] + BoardCoords.size(0)) - 1] -
                  BoardCoords[r3[i] - 1];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                x[i] = BoardCoords[(r2[i] + BoardCoords.size(0)) - 1] -
                  BoardCoords[r3[i] - 1];
              }
            }

            end = r.size(0) - 1;
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (r[b_i] && r1[b_i]) {
                loop_ub++;
              }
            }

            r4.set_size(loop_ub);
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (r[b_i] && r1[b_i]) {
                r4[loop_ub] = b_i + 1;
                loop_ub++;
              }
            }

            end = r.size(0) - 1;
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (r[b_i] && r1[b_i]) {
                loop_ub++;
              }
            }

            r5.set_size(loop_ub);
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (r[b_i] && r1[b_i]) {
                r5[loop_ub] = b_i + 1;
                loop_ub++;
              }
            }

            b_x.set_size(r4.size(0));
            loop_ub = r4.size(0);
            if (static_cast<int>(r4.size(0) < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                b_x[i] = BoardCoords[((r4[i] + BoardCoords.size(0)) +
                                      BoardCoords.size(0) * BoardCoords.size(1))
                  - 1] - BoardCoords[(r5[i] + BoardCoords.size(0) *
                                      BoardCoords.size(1)) - 1];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                b_x[i] = BoardCoords[((r4[i] + BoardCoords.size(0)) +
                                      BoardCoords.size(0) * BoardCoords.size(1))
                  - 1] - BoardCoords[(r5[i] + BoardCoords.size(0) *
                                      BoardCoords.size(1)) - 1];
              }
            }

            if (std::abs(combineVectorElements(x) / static_cast<double>(x.size(0)))
                > std::abs(combineVectorElements(b_x) / static_cast<double>
                           (b_x.size(0)))) {
              coordsToUse[0] = 1.0;
              coordsToUse[1] = 2.0;
            } else {
              coordsToUse[0] = 2.0;
              coordsToUse[1] = 1.0;
            }
          }

          //
          // Arguments    : const ::coder::array<double, 2U> &idx
          //                double coordsToUse[2]
          // Return Type  : void
          //
          void Checkerboard::b_findIndependentVar(const ::coder::array<double,
            2U> &idx, double coordsToUse[2]) const
          {
            array<double, 1U> b_x;
            array<double, 1U> x;
            array<int, 1U> r2;
            array<int, 1U> r3;
            array<int, 1U> r4;
            array<int, 1U> r5;
            array<boolean_T, 1U> r;
            array<boolean_T, 1U> r1;
            int end;
            int idx_tmp;
            int loop_ub;
            idx_tmp = static_cast<int>(idx[0]);
            r.set_size(BoardIdx.size(0));
            loop_ub = BoardIdx.size(0);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                r[i] = (BoardIdx[i + BoardIdx.size(0) * (idx_tmp - 1)] > 0.0);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                r[i] = (BoardIdx[i + BoardIdx.size(0) * (idx_tmp - 1)] > 0.0);
              }
            }

            idx_tmp = static_cast<int>(idx[1]);
            r1.set_size(BoardIdx.size(0));
            loop_ub = BoardIdx.size(0);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                r1[i] = (BoardIdx[i + BoardIdx.size(0) * (idx_tmp - 1)] > 0.0);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                r1[i] = (BoardIdx[i + BoardIdx.size(0) * (idx_tmp - 1)] > 0.0);
              }
            }

            end = r.size(0) - 1;
            idx_tmp = 0;
            for (loop_ub = 0; loop_ub <= end; loop_ub++) {
              if (r[loop_ub] && r1[loop_ub]) {
                idx_tmp++;
              }
            }

            r2.set_size(idx_tmp);
            idx_tmp = 0;
            for (loop_ub = 0; loop_ub <= end; loop_ub++) {
              if (r[loop_ub] && r1[loop_ub]) {
                r2[idx_tmp] = loop_ub + 1;
                idx_tmp++;
              }
            }

            end = r.size(0) - 1;
            idx_tmp = 0;
            for (loop_ub = 0; loop_ub <= end; loop_ub++) {
              if (r[loop_ub] && r1[loop_ub]) {
                idx_tmp++;
              }
            }

            r3.set_size(idx_tmp);
            idx_tmp = 0;
            for (loop_ub = 0; loop_ub <= end; loop_ub++) {
              if (r[loop_ub] && r1[loop_ub]) {
                r3[idx_tmp] = loop_ub + 1;
                idx_tmp++;
              }
            }

            idx_tmp = static_cast<int>(idx[1]);
            end = static_cast<int>(idx[0]);
            x.set_size(r2.size(0));
            loop_ub = r2.size(0);
            if (static_cast<int>(r2.size(0) < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                x[i] = BoardCoords[(r2[i] + BoardCoords.size(0) * (idx_tmp - 1))
                  - 1] - BoardCoords[(r3[i] + BoardCoords.size(0) * (end - 1)) -
                  1];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                x[i] = BoardCoords[(r2[i] + BoardCoords.size(0) * (idx_tmp - 1))
                  - 1] - BoardCoords[(r3[i] + BoardCoords.size(0) * (end - 1)) -
                  1];
              }
            }

            end = r.size(0) - 1;
            idx_tmp = 0;
            for (loop_ub = 0; loop_ub <= end; loop_ub++) {
              if (r[loop_ub] && r1[loop_ub]) {
                idx_tmp++;
              }
            }

            r4.set_size(idx_tmp);
            idx_tmp = 0;
            for (loop_ub = 0; loop_ub <= end; loop_ub++) {
              if (r[loop_ub] && r1[loop_ub]) {
                r4[idx_tmp] = loop_ub + 1;
                idx_tmp++;
              }
            }

            end = r.size(0) - 1;
            idx_tmp = 0;
            for (loop_ub = 0; loop_ub <= end; loop_ub++) {
              if (r[loop_ub] && r1[loop_ub]) {
                idx_tmp++;
              }
            }

            r5.set_size(idx_tmp);
            idx_tmp = 0;
            for (loop_ub = 0; loop_ub <= end; loop_ub++) {
              if (r[loop_ub] && r1[loop_ub]) {
                r5[idx_tmp] = loop_ub + 1;
                idx_tmp++;
              }
            }

            idx_tmp = static_cast<int>(idx[1]);
            end = static_cast<int>(idx[0]);
            b_x.set_size(r4.size(0));
            loop_ub = r4.size(0);
            if (static_cast<int>(r4.size(0) < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                b_x[i] = BoardCoords[((r4[i] + BoardCoords.size(0) * (idx_tmp -
                  1)) + BoardCoords.size(0) * BoardCoords.size(1)) - 1] -
                  BoardCoords[((r5[i] + BoardCoords.size(0) * (end - 1)) +
                               BoardCoords.size(0) * BoardCoords.size(1)) - 1];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                b_x[i] = BoardCoords[((r4[i] + BoardCoords.size(0) * (idx_tmp -
                  1)) + BoardCoords.size(0) * BoardCoords.size(1)) - 1] -
                  BoardCoords[((r5[i] + BoardCoords.size(0) * (end - 1)) +
                               BoardCoords.size(0) * BoardCoords.size(1)) - 1];
              }
            }

            if (std::abs(combineVectorElements(x) / static_cast<double>(x.size(0)))
                > std::abs(combineVectorElements(b_x) / static_cast<double>
                           (b_x.size(0)))) {
              coordsToUse[0] = 1.0;
              coordsToUse[1] = 2.0;
            } else {
              coordsToUse[0] = 2.0;
              coordsToUse[1] = 1.0;
            }
          }

          //
          // Arguments    : const ::coder::array<double, 2U> &idx
          //                ::coder::array<double, 2U> &newIndices
          // Return Type  : void
          //
          void Checkerboard::b_fitPolynomialIndices(const ::coder::array<double,
            2U> &idx, ::coder::array<double, 2U> &newIndices) const
          {
            array<double, 2U> b_this;
            array<double, 2U> removedIdx;
            array<double, 2U> validIdx;
            array<int, 2U> r;
            array<int, 1U> currCurve_tmp;
            double currCurve_data[5];
            double coordsToUse[2];
            double coordDist;
            double moveDistMultiplier;
            double refCoordValue;
            int b_i;
            int loop_ub;
            int n;
            b_findIndependentVar(idx, coordsToUse);
            newIndices.set_size(1, BoardCoords.size(0));
            loop_ub = BoardCoords.size(0);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                newIndices[i] = 0.0;
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                newIndices[i] = 0.0;
              }
            }

            removedIdx.set_size(1, 0);
            b_i = BoardCoords.size(0);
            for (int j{0}; j < b_i; j++) {
              int i1;
              validIdx.set_size(1, BoardCoords.size(1));
              loop_ub = BoardCoords.size(1);
              for (i1 = 0; i1 < loop_ub; i1++) {
                validIdx[i1] = BoardCoords[(j + BoardCoords.size(0) * i1) +
                  BoardCoords.size(0) * BoardCoords.size(1) * (static_cast<int>
                  (coordsToUse[0]) - 1)];
              }

              eml_find(validIdx, r);
              validIdx.set_size(1, r.size(1));
              loop_ub = r.size(1);
              for (i1 = 0; i1 < loop_ub; i1++) {
                validIdx[i1] = r[i1];
              }

              if (validIdx.size(1) >= 2) {
                double currCoord;
                double currRad;
                int currCurve_size[2];
                boolean_T exitg1;
                findSearchParams(idx, validIdx, static_cast<double>(j) + 1.0,
                                 coordsToUse, &coordDist, &moveDistMultiplier,
                                 &refCoordValue);
                n = 0;
                i1 = validIdx.size(1);
                for (loop_ub = 0; loop_ub < i1; loop_ub++) {
                  if (static_cast<int>(validIdx[loop_ub]) != 0) {
                    n++;
                  }
                }

                currCurve_tmp.set_size(validIdx.size(1));
                loop_ub = validIdx.size(1);
                for (i1 = 0; i1 < loop_ub; i1++) {
                  currCurve_tmp[i1] = static_cast<int>(validIdx[i1]);
                }

                validIdx.set_size(1, currCurve_tmp.size(0));
                loop_ub = currCurve_tmp.size(0);
                for (i1 = 0; i1 < loop_ub; i1++) {
                  validIdx[i1] = BoardCoords[(j + BoardCoords.size(0) *
                    (currCurve_tmp[i1] - 1)) + BoardCoords.size(0) *
                    BoardCoords.size(1) * (static_cast<int>(coordsToUse[0]) - 1)];
                }

                b_this.set_size(1, currCurve_tmp.size(0));
                loop_ub = currCurve_tmp.size(0);
                for (i1 = 0; i1 < loop_ub; i1++) {
                  b_this[i1] = BoardCoords[(j + BoardCoords.size(0) *
                    (currCurve_tmp[i1] - 1)) + BoardCoords.size(0) *
                    BoardCoords.size(1) * (static_cast<int>(coordsToUse[1]) - 1)];
                }

                if (n > 5) {
                  i1 = 4;
                } else {
                  i1 = 2;
                }

                polyfit(validIdx, b_this, static_cast<double>(i1),
                        currCurve_data, currCurve_size);
                currRad = coordDist / 4.0;
                refCoordValue = BoardCoords[(j + BoardCoords.size(0) * (
                  static_cast<int>(refCoordValue) - 1)) + BoardCoords.size(0) *
                  BoardCoords.size(1) * (static_cast<int>(coordsToUse[0]) - 1)];
                currCoord = currRad + refCoordValue;
                exitg1 = false;
                while ((!exitg1) && (std::abs(currCoord - refCoordValue) <
                                     moveDistMultiplier * 1.5 * std::abs
                                     (coordDist))) {
                  double currPt[2];
                  boolean_T exitg2;
                  boolean_T p;
                  p = true;
                  loop_ub = 0;
                  exitg2 = false;
                  while ((!exitg2) && (loop_ub < 2)) {
                    if (!(coordsToUse[loop_ub] == static_cast<double>(loop_ub) +
                          1.0)) {
                      p = false;
                      exitg2 = true;
                    } else {
                      loop_ub++;
                    }
                  }

                  if (p) {
                    double y;
                    y = currCurve_data[0];
                    i1 = currCurve_size[1];
                    for (loop_ub = 0; loop_ub <= i1 - 2; loop_ub++) {
                      y = currCoord * y + currCurve_data[loop_ub + 1];
                    }

                    currPt[0] = currCoord;
                    currPt[1] = y;
                  } else {
                    double y;
                    y = currCurve_data[0];
                    i1 = currCurve_size[1];
                    for (loop_ub = 0; loop_ub <= i1 - 2; loop_ub++) {
                      y = currCoord * y + currCurve_data[loop_ub + 1];
                    }

                    currPt[0] = y;
                    currPt[1] = currCoord;
                  }

                  findClosestOnCurve(currPt, std::abs(currRad), currCurve_data,
                                     currCurve_size, coordsToUse, removedIdx,
                                     validIdx);
                  if (validIdx.size(1) != 0) {
                    newIndices[j] = validIdx[0];
                    i1 = removedIdx.size(1);
                    loop_ub = validIdx.size(1);
                    removedIdx.set_size(removedIdx.size(0), removedIdx.size(1) +
                                        validIdx.size(1));
                    for (n = 0; n < loop_ub; n++) {
                      removedIdx[i1 + n] = validIdx[n];
                    }

                    exitg1 = true;
                  } else {
                    currCoord += currRad;
                  }
                }
              }
            }

            n = 0;
            b_i = newIndices.size(1);
            for (loop_ub = 0; loop_ub < b_i; loop_ub++) {
              if (newIndices[loop_ub] != 0.0) {
                n++;
              }
            }

            if (n < 4) {
              n = newIndices.size(1);
              if (static_cast<int>(newIndices.size(1) < 3200)) {
                for (int i{0}; i < n; i++) {
                  if (newIndices[i] > 0.0) {
                    newIndices[i] = 0.0;
                  }
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < n; i++) {
                  if (newIndices[i] > 0.0) {
                    newIndices[i] = 0.0;
                  }
                }
              }
            }
          }

          //
          // Arguments    : ::coder::array<double, 2U> &newIndices
          // Return Type  : void
          //
          void Checkerboard::b_fitPolynomialIndices(::coder::array<double, 2U>
            &newIndices) const
          {
            array<double, 2U> b_index;
            array<double, 2U> b_this;
            array<double, 2U> removedIdx;
            array<int, 2U> validIdx;
            array<int, 1U> currCurve_tmp;
            double currCurve_data[5];
            double coordsToUse[2];
            int i1;
            int loop_ub;
            b_findIndependentVar(coordsToUse);
            newIndices.set_size(1, BoardCoords.size(0));
            loop_ub = BoardCoords.size(0);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                newIndices[i] = 0.0;
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                newIndices[i] = 0.0;
              }
            }

            removedIdx.set_size(1, 0);
            i1 = BoardCoords.size(0);
            for (int j{0}; j < i1; j++) {
              int i2;
              b_index.set_size(1, BoardCoords.size(1));
              loop_ub = BoardCoords.size(1);
              for (i2 = 0; i2 < loop_ub; i2++) {
                b_index[i2] = BoardCoords[(j + BoardCoords.size(0) * i2) +
                  BoardCoords.size(0) * BoardCoords.size(1) * (static_cast<int>
                  (coordsToUse[0]) - 1)];
              }

              eml_find(b_index, validIdx);
              if (validIdx.size(1) >= 2) {
                double coordDist;
                double currCoord;
                double currRad;
                double refCoordValue;
                int currCurve_size[2];
                int coordDist_tmp;
                int i3;
                int n;
                boolean_T exitg1;
                coordDist_tmp = validIdx[0];
                coordDist = (BoardCoords[(j + BoardCoords.size(0) * (validIdx[0]
                  - 1)) + BoardCoords.size(0) * BoardCoords.size(1) * (
                  static_cast<int>(coordsToUse[0]) - 1)] - BoardCoords[(j +
                  BoardCoords.size(0) * (validIdx[1] - 1)) + BoardCoords.size(0)
                             * BoardCoords.size(1) * (static_cast<int>
                  (coordsToUse[0]) - 1)]) / (static_cast<double>(validIdx[1]) -
                  static_cast<double>(coordDist_tmp));
                n = 0;
                i2 = validIdx.size(1);
                currCurve_tmp.set_size(validIdx.size(1));
                for (loop_ub = 0; loop_ub < i2; loop_ub++) {
                  i3 = validIdx[loop_ub];
                  if (i3 != 0) {
                    n++;
                  }

                  currCurve_tmp[loop_ub] = i3;
                }

                b_index.set_size(1, currCurve_tmp.size(0));
                loop_ub = currCurve_tmp.size(0);
                for (i2 = 0; i2 < loop_ub; i2++) {
                  b_index[i2] = BoardCoords[(j + BoardCoords.size(0) *
                    (currCurve_tmp[i2] - 1)) + BoardCoords.size(0) *
                    BoardCoords.size(1) * (static_cast<int>(coordsToUse[0]) - 1)];
                }

                b_this.set_size(1, currCurve_tmp.size(0));
                loop_ub = currCurve_tmp.size(0);
                for (i2 = 0; i2 < loop_ub; i2++) {
                  b_this[i2] = BoardCoords[(j + BoardCoords.size(0) *
                    (currCurve_tmp[i2] - 1)) + BoardCoords.size(0) *
                    BoardCoords.size(1) * (static_cast<int>(coordsToUse[1]) - 1)];
                }

                if (n > 5) {
                  i2 = 4;
                } else {
                  i2 = 2;
                }

                polyfit(b_index, b_this, static_cast<double>(i2), currCurve_data,
                        currCurve_size);
                currRad = coordDist / 4.0;
                refCoordValue = BoardCoords[(j + BoardCoords.size(0) *
                  (validIdx[0] - 1)) + BoardCoords.size(0) * BoardCoords.size(1)
                  * (static_cast<int>(coordsToUse[0]) - 1)];
                currCoord = currRad + refCoordValue;
                exitg1 = false;
                while ((!exitg1) && (std::abs(currCoord - refCoordValue) <
                                     static_cast<double>(coordDist_tmp) * 1.5 *
                                     std::abs(coordDist))) {
                  double currPt[2];
                  boolean_T exitg2;
                  boolean_T p;
                  p = true;
                  loop_ub = 0;
                  exitg2 = false;
                  while ((!exitg2) && (loop_ub < 2)) {
                    if (!(coordsToUse[loop_ub] == static_cast<double>(loop_ub) +
                          1.0)) {
                      p = false;
                      exitg2 = true;
                    } else {
                      loop_ub++;
                    }
                  }

                  if (p) {
                    double y;
                    y = currCurve_data[0];
                    i2 = currCurve_size[1];
                    for (loop_ub = 0; loop_ub <= i2 - 2; loop_ub++) {
                      y = currCoord * y + currCurve_data[loop_ub + 1];
                    }

                    currPt[0] = currCoord;
                    currPt[1] = y;
                  } else {
                    double y;
                    y = currCurve_data[0];
                    i2 = currCurve_size[1];
                    for (loop_ub = 0; loop_ub <= i2 - 2; loop_ub++) {
                      y = currCoord * y + currCurve_data[loop_ub + 1];
                    }

                    currPt[0] = y;
                    currPt[1] = currCoord;
                  }

                  findClosestOnCurve(currPt, std::abs(currRad), currCurve_data,
                                     currCurve_size, coordsToUse, removedIdx,
                                     b_index);
                  if (b_index.size(1) != 0) {
                    newIndices[j] = b_index[0];
                    i2 = removedIdx.size(1);
                    loop_ub = b_index.size(1);
                    removedIdx.set_size(removedIdx.size(0), removedIdx.size(1) +
                                        b_index.size(1));
                    for (i3 = 0; i3 < loop_ub; i3++) {
                      removedIdx[i2 + i3] = b_index[i3];
                    }

                    exitg1 = true;
                  } else {
                    currCoord += currRad;
                  }
                }
              }
            }
          }

          //
          // Arguments    : ::coder::array<double, 2U> &newIndices
          // Return Type  : void
          //
          void Checkerboard::c_fitPolynomialIndices(::coder::array<double, 2U>
            &newIndices) const
          {
            array<double, 2U> b_index;
            array<double, 2U> removedIdx;
            array<double, 1U> b_this;
            array<double, 1U> c_this;
            array<int, 1U> validIdx;
            double currCurve_data[5];
            double coordsToUse[2];
            int b_i;
            int k;
            int loop_ub;
            findIndependentVar(coordsToUse);
            newIndices.set_size(1, BoardCoords.size(1));
            loop_ub = BoardCoords.size(1);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                newIndices[i] = 0.0;
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                newIndices[i] = 0.0;
              }
            }

            removedIdx.set_size(1, 0);
            b_i = BoardCoords.size(1);
            for (int j{0}; j < b_i; j++) {
              int i1;
              b_this.set_size(BoardCoords.size(0));
              loop_ub = BoardCoords.size(0);
              for (i1 = 0; i1 < loop_ub; i1++) {
                b_this[i1] = BoardCoords[(i1 + BoardCoords.size(0) * j) +
                  BoardCoords.size(0) * BoardCoords.size(1) * (static_cast<int>
                  (coordsToUse[0]) - 1)];
              }

              eml_find(b_this, validIdx);
              if (validIdx.size(0) >= 2) {
                double coordDist;
                double currCoord;
                double currRad;
                double refCoordValue;
                int currCurve_size[2];
                boolean_T exitg1;
                coordDist = (BoardCoords[((validIdx[0] + BoardCoords.size(0) * j)
                  + BoardCoords.size(0) * BoardCoords.size(1) * (static_cast<int>
                  (coordsToUse[0]) - 1)) - 1] - BoardCoords[((validIdx[1] +
                  BoardCoords.size(0) * j) + BoardCoords.size(0) *
                  BoardCoords.size(1) * (static_cast<int>(coordsToUse[0]) - 1))
                             - 1]) / (static_cast<double>(validIdx[1]) -
                                      static_cast<double>(validIdx[0]));
                loop_ub = 0;
                i1 = validIdx.size(0);
                b_this.set_size(validIdx.size(0));
                c_this.set_size(validIdx.size(0));
                for (k = 0; k < i1; k++) {
                  if (validIdx[k] != 0) {
                    loop_ub++;
                  }

                  b_this[k] = BoardCoords[((validIdx[k] + BoardCoords.size(0) *
                    j) + BoardCoords.size(0) * BoardCoords.size(1) * (
                    static_cast<int>(coordsToUse[0]) - 1)) - 1];
                  c_this[k] = BoardCoords[((validIdx[k] + BoardCoords.size(0) *
                    j) + BoardCoords.size(0) * BoardCoords.size(1) * (
                    static_cast<int>(coordsToUse[1]) - 1)) - 1];
                }

                if (loop_ub > 5) {
                  i1 = 4;
                } else {
                  i1 = 2;
                }

                polyfit(b_this, c_this, static_cast<double>(i1), currCurve_data,
                        currCurve_size);
                currRad = coordDist / 4.0;
                refCoordValue = BoardCoords[((validIdx[0] + BoardCoords.size(0) *
                  j) + BoardCoords.size(0) * BoardCoords.size(1) * (static_cast<
                  int>(coordsToUse[0]) - 1)) - 1];
                currCoord = currRad + refCoordValue;
                exitg1 = false;
                while ((!exitg1) && (std::abs(currCoord - refCoordValue) <
                                     static_cast<double>(validIdx[0]) * 1.5 *
                                     std::abs(coordDist))) {
                  double currPt[2];
                  boolean_T exitg2;
                  boolean_T p;
                  p = true;
                  k = 0;
                  exitg2 = false;
                  while ((!exitg2) && (k < 2)) {
                    if (!(coordsToUse[k] == static_cast<double>(k) + 1.0)) {
                      p = false;
                      exitg2 = true;
                    } else {
                      k++;
                    }
                  }

                  if (p) {
                    double y;
                    y = currCurve_data[0];
                    i1 = currCurve_size[1];
                    for (k = 0; k <= i1 - 2; k++) {
                      y = currCoord * y + currCurve_data[k + 1];
                    }

                    currPt[0] = currCoord;
                    currPt[1] = y;
                  } else {
                    double y;
                    y = currCurve_data[0];
                    i1 = currCurve_size[1];
                    for (k = 0; k <= i1 - 2; k++) {
                      y = currCoord * y + currCurve_data[k + 1];
                    }

                    currPt[0] = y;
                    currPt[1] = currCoord;
                  }

                  findClosestOnCurve(currPt, std::abs(currRad), currCurve_data,
                                     currCurve_size, coordsToUse, removedIdx,
                                     b_index);
                  if (b_index.size(1) != 0) {
                    newIndices[j] = b_index[0];
                    i1 = removedIdx.size(1);
                    loop_ub = b_index.size(1);
                    removedIdx.set_size(removedIdx.size(0), removedIdx.size(1) +
                                        b_index.size(1));
                    for (k = 0; k < loop_ub; k++) {
                      removedIdx[i1 + k] = b_index[k];
                    }

                    exitg1 = true;
                  } else {
                    currCoord += currRad;
                  }
                }
              }
            }

            loop_ub = 0;
            b_i = newIndices.size(1);
            for (k = 0; k < b_i; k++) {
              if (newIndices[k] != 0.0) {
                loop_ub++;
              }
            }

            if (loop_ub < 4) {
              loop_ub = newIndices.size(1);
              if (static_cast<int>(newIndices.size(1) < 3200)) {
                for (int i{0}; i < loop_ub; i++) {
                  if (newIndices[i] > 0.0) {
                    newIndices[i] = 0.0;
                  }
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < loop_ub; i++) {
                  if (newIndices[i] > 0.0) {
                    newIndices[i] = 0.0;
                  }
                }
              }
            }
          }

          //
          // Arguments    : void
          // Return Type  : float
          //
          float Checkerboard::computeInitialEnergy() const
          {
            array<float, 2U> col1;
            array<float, 2U> col2;
            array<float, 2U> row3;
            array<boolean_T, 1U> x;
            float e;
            int loop_ub;
            boolean_T exitg1;
            boolean_T y;
            x.set_size(BoardIdx.size(0) * BoardIdx.size(1));
            loop_ub = BoardIdx.size(0) * BoardIdx.size(1);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                x[i] = (BoardIdx[i] < 0.0);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                x[i] = (BoardIdx[i] < 0.0);
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

            if (y) {
              e = rtInfF;
            } else {
              float b_x[3];
              float c_x[3];
              float boardSize;
              float y_idx_0;
              float y_idx_1;
              float y_idx_2;
              col1.set_size(3, Points.size(1));
              loop_ub = Points.size(1);
              if (static_cast<int>(loop_ub * 3 < 3200)) {
                for (int i{0}; i < loop_ub; i++) {
                  col1[3 * i] = Points[(static_cast<int>(BoardIdx[0]) +
                                        Points.size(0) * i) - 1];
                  col1[3 * i + 1] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0)]) + Points.size(0) * i) - 1];
                  col1[3 * i + 2] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0) * 2]) + Points.size(0) * i) - 1];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < loop_ub; i++) {
                  col1[3 * i] = Points[(static_cast<int>(BoardIdx[0]) +
                                        Points.size(0) * i) - 1];
                  col1[3 * i + 1] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0)]) + Points.size(0) * i) - 1];
                  col1[3 * i + 2] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0) * 2]) + Points.size(0) * i) - 1];
                }
              }

              col2.set_size(3, Points.size(1));
              loop_ub = Points.size(1);
              if (static_cast<int>(loop_ub * 3 < 3200)) {
                for (int i{0}; i < loop_ub; i++) {
                  col2[3 * i] = Points[(static_cast<int>(BoardIdx[1]) +
                                        Points.size(0) * i) - 1];
                  col2[3 * i + 1] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0) + 1]) + Points.size(0) * i) - 1];
                  col2[3 * i + 2] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0) * 2 + 1]) + Points.size(0) * i) -
                    1];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < loop_ub; i++) {
                  col2[3 * i] = Points[(static_cast<int>(BoardIdx[1]) +
                                        Points.size(0) * i) - 1];
                  col2[3 * i + 1] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0) + 1]) + Points.size(0) * i) - 1];
                  col2[3 * i + 2] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0) * 2 + 1]) + Points.size(0) * i) -
                    1];
                }
              }

              row3.set_size(3, Points.size(1));
              loop_ub = Points.size(1);
              if (static_cast<int>(loop_ub * 3 < 3200)) {
                for (int i{0}; i < loop_ub; i++) {
                  row3[3 * i] = Points[(static_cast<int>(BoardIdx[2]) +
                                        Points.size(0) * i) - 1];
                  row3[3 * i + 1] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0) + 2]) + Points.size(0) * i) - 1];
                  row3[3 * i + 2] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0) * 2 + 2]) + Points.size(0) * i) -
                    1];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < loop_ub; i++) {
                  row3[3 * i] = Points[(static_cast<int>(BoardIdx[2]) +
                                        Points.size(0) * i) - 1];
                  row3[3 * i + 1] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0) + 2]) + Points.size(0) * i) - 1];
                  row3[3 * i + 2] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0) * 2 + 2]) + Points.size(0) * i) -
                    1];
                }
              }

              if (col1.size(1) == 1) {
                loop_ub = row3.size(1);
              } else {
                loop_ub = col1.size(1);
              }

              if ((col1.size(1) == row3.size(1)) && (loop_ub == col2.size(1))) {
                loop_ub = 3 * col1.size(1);
                col2.set_size(3, col1.size(1));
                if (static_cast<int>(loop_ub < 3200)) {
                  for (int i{0}; i < loop_ub; i++) {
                    col2[i] = (col1[i] + row3[i]) - 2.0F * col2[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (int i = 0; i < loop_ub; i++) {
                    col2[i] = (col1[i] + row3[i]) - 2.0F * col2[i];
                  }
                }
              } else {
                c_binary_expand_op(col2, col1, row3);
              }

              if (col1.size(1) == row3.size(1)) {
                loop_ub = 3 * col1.size(1);
                col1.set_size(3, col1.size(1));
                if (static_cast<int>(loop_ub < 3200)) {
                  for (int i{0}; i < loop_ub; i++) {
                    col1[i] = col1[i] - row3[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (int i = 0; i < loop_ub; i++) {
                    col1[i] = col1[i] - row3[i];
                  }
                }
              } else {
                b_minus(col1, row3);
              }

              b_x[0] = rt_hypotf_snf(col2[0], col2[3]);
              y_idx_0 = rt_hypotf_snf(col1[0], col1[3]);
              b_x[1] = rt_hypotf_snf(col2[1], col2[4]);
              y_idx_1 = rt_hypotf_snf(col1[1], col1[4]);
              b_x[2] = rt_hypotf_snf(col2[2], col2[5]);
              y_idx_2 = rt_hypotf_snf(col1[2], col1[5]);
              col1.set_size(3, Points.size(1));
              loop_ub = Points.size(1);
              if (static_cast<int>(loop_ub * 3 < 3200)) {
                for (int i{0}; i < loop_ub; i++) {
                  col1[3 * i] = Points[(static_cast<int>(BoardIdx[0]) +
                                        Points.size(0) * i) - 1];
                  col1[3 * i + 1] = Points[(static_cast<int>(BoardIdx[1]) +
                    Points.size(0) * i) - 1];
                  col1[3 * i + 2] = Points[(static_cast<int>(BoardIdx[2]) +
                    Points.size(0) * i) - 1];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < loop_ub; i++) {
                  col1[3 * i] = Points[(static_cast<int>(BoardIdx[0]) +
                                        Points.size(0) * i) - 1];
                  col1[3 * i + 1] = Points[(static_cast<int>(BoardIdx[1]) +
                    Points.size(0) * i) - 1];
                  col1[3 * i + 2] = Points[(static_cast<int>(BoardIdx[2]) +
                    Points.size(0) * i) - 1];
                }
              }

              col2.set_size(3, Points.size(1));
              loop_ub = Points.size(1);
              if (static_cast<int>(loop_ub * 3 < 3200)) {
                for (int i{0}; i < loop_ub; i++) {
                  col2[3 * i] = Points[(static_cast<int>(BoardIdx[BoardIdx.size
                    (0)]) + Points.size(0) * i) - 1];
                  col2[3 * i + 1] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0) + 1]) + Points.size(0) * i) - 1];
                  col2[3 * i + 2] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0) + 2]) + Points.size(0) * i) - 1];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < loop_ub; i++) {
                  col2[3 * i] = Points[(static_cast<int>(BoardIdx[BoardIdx.size
                    (0)]) + Points.size(0) * i) - 1];
                  col2[3 * i + 1] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0) + 1]) + Points.size(0) * i) - 1];
                  col2[3 * i + 2] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0) + 2]) + Points.size(0) * i) - 1];
                }
              }

              row3.set_size(3, Points.size(1));
              loop_ub = Points.size(1);
              if (static_cast<int>(loop_ub * 3 < 3200)) {
                for (int i{0}; i < loop_ub; i++) {
                  row3[3 * i] = Points[(static_cast<int>(BoardIdx[BoardIdx.size
                    (0) * 2]) + Points.size(0) * i) - 1];
                  row3[3 * i + 1] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0) * 2 + 1]) + Points.size(0) * i) -
                    1];
                  row3[3 * i + 2] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0) * 2 + 2]) + Points.size(0) * i) -
                    1];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < loop_ub; i++) {
                  row3[3 * i] = Points[(static_cast<int>(BoardIdx[BoardIdx.size
                    (0) * 2]) + Points.size(0) * i) - 1];
                  row3[3 * i + 1] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0) * 2 + 1]) + Points.size(0) * i) -
                    1];
                  row3[3 * i + 2] = Points[(static_cast<int>
                    (BoardIdx[BoardIdx.size(0) * 2 + 2]) + Points.size(0) * i) -
                    1];
                }
              }

              if (col1.size(1) == 1) {
                loop_ub = row3.size(1);
              } else {
                loop_ub = col1.size(1);
              }

              if ((col1.size(1) == row3.size(1)) && (loop_ub == col2.size(1))) {
                loop_ub = 3 * col1.size(1);
                col2.set_size(3, col1.size(1));
                if (static_cast<int>(loop_ub < 3200)) {
                  for (int i{0}; i < loop_ub; i++) {
                    col2[i] = (col1[i] + row3[i]) - 2.0F * col2[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (int i = 0; i < loop_ub; i++) {
                    col2[i] = (col1[i] + row3[i]) - 2.0F * col2[i];
                  }
                }
              } else {
                c_binary_expand_op(col2, col1, row3);
              }

              if (col1.size(1) == row3.size(1)) {
                loop_ub = 3 * col1.size(1);
                col1.set_size(3, col1.size(1));
                if (static_cast<int>(loop_ub < 3200)) {
                  for (int i{0}; i < loop_ub; i++) {
                    col1[i] = col1[i] - row3[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (int i = 0; i < loop_ub; i++) {
                    col1[i] = col1[i] - row3[i];
                  }
                }
              } else {
                b_minus(col1, row3);
              }

              boardSize = static_cast<float>(BoardIdx.size(0) * BoardIdx.size(1));
              b_x[0] /= y_idx_0;
              c_x[0] = rt_hypotf_snf(col2[0], col2[3]) / rt_hypotf_snf(col1[0],
                col1[3]);
              b_x[1] /= y_idx_1;
              c_x[1] = rt_hypotf_snf(col2[1], col2[4]) / rt_hypotf_snf(col1[1],
                col1[4]);
              b_x[2] /= y_idx_2;
              c_x[2] = rt_hypotf_snf(col2[2], col2[5]) / rt_hypotf_snf(col1[2],
                col1[5]);
              e = boardSize * std::fmax(std::fmax(0.0F, ::coder::internal::
                maximum(b_x)), ::coder::internal::maximum(c_x)) - boardSize;
            }

            return e;
          }

          //
          // Arguments    : float oldEnergy
          // Return Type  : float
          //
          float Checkerboard::computeNewEnergyHorizontal(float oldEnergy) const
          {
            array<double, 3U> c_this;
            array<double, 3U> denom;
            array<double, 3U> num;
            array<double, 2U> b_denom;
            array<double, 2U> b_num;
            array<double, 2U> validNewColIdx;
            array<double, 1U> b_r;
            array<double, 1U> c_r;
            array<int, 1U> r;
            array<int, 1U> r1;
            array<int, 1U> r2;
            array<int, 1U> r3;
            array<int, 1U> r4;
            array<boolean_T, 2U> b_this;
            array<boolean_T, 1U> validIdx;
            float newEnergy;
            int end;
            int i2;
            int k;
            int loop_ub;
            int nx;
            boolean_T exitg1;
            boolean_T y;
            validIdx.set_size(BoardIdx.size(0));
            loop_ub = BoardIdx.size(0);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                validIdx[i] = ((BoardIdx[i] > 0.0) && (BoardIdx[i +
                  BoardIdx.size(0)] > 0.0) && (BoardIdx[i + BoardIdx.size(0) * 2]
                  > 0.0));
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                validIdx[i] = ((BoardIdx[i] > 0.0) && (BoardIdx[i +
                  BoardIdx.size(0)] > 0.0) && (BoardIdx[i + BoardIdx.size(0) * 2]
                  > 0.0));
              }
            }

            newEnergy = 0.0F;
            y = false;
            nx = 1;
            exitg1 = false;
            while ((!exitg1) && (nx <= validIdx.size(0))) {
              if (validIdx[nx - 1]) {
                y = true;
                exitg1 = true;
              } else {
                nx++;
              }
            }

            if (y) {
              end = validIdx.size(0) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r.set_size(nx);
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  r[nx] = b_i + 1;
                  nx++;
                }
              }

              end = validIdx.size(0) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r1.set_size(nx);
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  r1[nx] = b_i + 1;
                  nx++;
                }
              }

              end = validIdx.size(0) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r2.set_size(nx);
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  r2[nx] = b_i + 1;
                  nx++;
                }
              }

              if (r.size(0) == r2.size(0)) {
                c_this.set_size(r.size(0), 1, BoardCoords.size(2));
                loop_ub = BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i2,k)

                for (int i = 0; i < loop_ub; i++) {
                  k = r.size(0);
                  for (i2 = 0; i2 < k; i2++) {
                    c_this[i2 + c_this.size(0) * i] = (BoardCoords[(r[i2] +
                      BoardCoords.size(0) * BoardCoords.size(1) * i) - 1] +
                      BoardCoords[((r1[i2] + BoardCoords.size(0) * 2) +
                                   BoardCoords.size(0) * BoardCoords.size(1) * i)
                      - 1]) - 2.0 * BoardCoords[((r2[i2] + BoardCoords.size(0))
                      + BoardCoords.size(0) * BoardCoords.size(1) * i) - 1];
                  }
                }

                b_squeeze(c_this, b_num);
              } else {
                binary_expand_op(b_num, this, r, r1, r2);
              }

              end = validIdx.size(0) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r3.set_size(nx);
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  r3[nx] = b_i + 1;
                  nx++;
                }
              }

              end = validIdx.size(0) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r4.set_size(nx);
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  r4[nx] = b_i + 1;
                  nx++;
                }
              }

              c_this.set_size(r3.size(0), 1, BoardCoords.size(2));
              loop_ub = BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i2,k)

              for (int i = 0; i < loop_ub; i++) {
                k = r3.size(0);
                for (i2 = 0; i2 < k; i2++) {
                  c_this[i2 + c_this.size(0) * i] = BoardCoords[(r3[i2] +
                    BoardCoords.size(0) * BoardCoords.size(1) * i) - 1] -
                    BoardCoords[((r4[i2] + BoardCoords.size(0) * 2) +
                                 BoardCoords.size(0) * BoardCoords.size(1) * i)
                    - 1];
                }
              }

              b_squeeze(c_this, b_denom);
              if (b_num.size(1) > 1) {
                b_r.set_size(b_num.size(0));
                nx = b_num.size(0);
                if (static_cast<int>(b_num.size(0) < 3200)) {
                  for (k = 0; k < nx; k++) {
                    b_r[k] = rt_hypotd_snf(b_num[k], b_num[k + b_num.size(0)]);
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (k = 0; k < nx; k++) {
                    b_r[k] = rt_hypotd_snf(b_num[k], b_num[k + b_num.size(0)]);
                  }
                }

                c_r.set_size(b_denom.size(0));
                nx = b_denom.size(0);
                if (static_cast<int>(b_denom.size(0) < 3200)) {
                  for (k = 0; k < nx; k++) {
                    c_r[k] = rt_hypotd_snf(b_denom[k], b_denom[k + b_denom.size
                      (0)]);
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (k = 0; k < nx; k++) {
                    c_r[k] = rt_hypotd_snf(b_denom[k], b_denom[k + b_denom.size
                      (0)]);
                  }
                }

                if (b_r.size(0) == c_r.size(0)) {
                  loop_ub = b_r.size(0);
                  if (static_cast<int>(b_r.size(0) < 3200)) {
                    for (int i{0}; i < loop_ub; i++) {
                      b_r[i] = b_r[i] / c_r[i];
                    }
                  } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                    for (int i = 0; i < loop_ub; i++) {
                      b_r[i] = b_r[i] / c_r[i];
                    }
                  }

                  newEnergy = std::fmax(oldEnergy, static_cast<float>(::coder::
                    internal::maximum(b_r)));
                } else {
                  newEnergy = binary_expand_op(oldEnergy, b_r, c_r);
                }
              } else {
                newEnergy = std::fmax(oldEnergy, static_cast<float>
                                      (rt_hypotd_snf(b_num[0], b_num[1]) /
                  rt_hypotd_snf(b_denom[0], b_denom[1])));
              }
            }

            b_this.set_size(1, BoardIdx.size(0));
            loop_ub = BoardIdx.size(0);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                b_this[i] = (BoardIdx[i] > 0.0);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                b_this[i] = (BoardIdx[i] > 0.0);
              }
            }

            Checkerboard::arrayFind(b_this, validNewColIdx);
            if (validNewColIdx.size(1) != 0) {
              int b_loop_ub;
              int i1;
              i1 = validNewColIdx.size(1);
              loop_ub = BoardCoords.size(2);
              b_loop_ub = BoardCoords.size(2);
              for (int b_i{0}; b_i < i1; b_i++) {
                double d;
                d = validNewColIdx[b_i];
                num.set_size(1, 1, BoardCoords.size(2));
                for (nx = 0; nx < loop_ub; nx++) {
                  num[nx] = (BoardCoords[(static_cast<int>(d) + BoardCoords.size
                              (0) * BoardCoords.size(1) * nx) - 1] +
                             BoardCoords[(static_cast<int>(d + 2.0) +
                              BoardCoords.size(0) * BoardCoords.size(1) * nx) -
                             1]) - 2.0 * BoardCoords[(static_cast<int>(d + 1.0)
                    + BoardCoords.size(0) * BoardCoords.size(1) * nx) - 1];
                }

                d = validNewColIdx[b_i];
                denom.set_size(1, 1, BoardCoords.size(2));
                for (nx = 0; nx < b_loop_ub; nx++) {
                  denom[nx] = BoardCoords[(static_cast<int>(d) +
                    BoardCoords.size(0) * BoardCoords.size(1) * nx) - 1] -
                    BoardCoords[(static_cast<int>(d + 2.0) + BoardCoords.size(0)
                                 * BoardCoords.size(1) * nx) - 1];
                }

                if (newEnergy != 0.0F) {
                  nx = num.size(2);
                  end = denom.size(2);
                  b_r = num.reshape(nx);
                  c_r = denom.reshape(end);
                  newEnergy = std::fmax(newEnergy, static_cast<float>(b_norm(b_r)
                    / b_norm(c_r)));
                } else {
                  nx = num.size(2);
                  end = denom.size(2);
                  b_r = num.reshape(nx);
                  c_r = denom.reshape(end);
                  newEnergy = std::fmax(oldEnergy, static_cast<float>(b_norm(b_r)
                    / b_norm(c_r)));
                }
              }
            }

            if (newEnergy != 0.0F) {
              newEnergy = newEnergy * static_cast<float>(BoardIdx.size(0) *
                BoardIdx.size(1)) - static_cast<float>(BoardIdx.size(0) *
                BoardIdx.size(1));
            } else {
              newEnergy = rtInfF;
            }

            return newEnergy;
          }

          //
          // Arguments    : const ::coder::array<double, 2U> &idx
          //                float oldEnergy
          // Return Type  : float
          //
          float Checkerboard::computeNewEnergyHorizontal(const ::coder::array<
            double, 2U> &idx, float oldEnergy) const
          {
            array<double, 3U> b_num;
            array<double, 3U> c_this;
            array<double, 3U> denom;
            array<double, 2U> b_denom;
            array<double, 2U> num;
            array<double, 2U> validNewColIdx;
            array<double, 1U> b_r;
            array<double, 1U> c_r;
            array<int, 1U> r;
            array<int, 1U> r1;
            array<int, 1U> r2;
            array<int, 1U> r3;
            array<int, 1U> r4;
            array<boolean_T, 2U> b_this;
            array<boolean_T, 1U> validIdx;
            float newEnergy;
            int b_idx;
            int b_idx_tmp;
            int c_idx;
            int d_idx;
            int i1;
            int idx_tmp;
            int k;
            int loop_ub;
            int nx;
            boolean_T exitg1;
            boolean_T y;
            nx = static_cast<int>(idx[0]);
            idx_tmp = static_cast<int>(idx[1]);
            b_idx_tmp = static_cast<int>(idx[2]);
            validIdx.set_size(BoardIdx.size(0));
            loop_ub = BoardIdx.size(0);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                validIdx[i] = ((BoardIdx[i + BoardIdx.size(0) * (nx - 1)] > 0.0)
                               && (BoardIdx[i + BoardIdx.size(0) * (idx_tmp - 1)]
                                   > 0.0) && (BoardIdx[i + BoardIdx.size(0) *
                  (b_idx_tmp - 1)] > 0.0));
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                validIdx[i] = ((BoardIdx[i + BoardIdx.size(0) * (nx - 1)] > 0.0)
                               && (BoardIdx[i + BoardIdx.size(0) * (idx_tmp - 1)]
                                   > 0.0) && (BoardIdx[i + BoardIdx.size(0) *
                  (b_idx_tmp - 1)] > 0.0));
              }
            }

            newEnergy = 0.0F;
            y = false;
            nx = 1;
            exitg1 = false;
            while ((!exitg1) && (nx <= validIdx.size(0))) {
              if (validIdx[nx - 1]) {
                y = true;
                exitg1 = true;
              } else {
                nx++;
              }
            }

            if (y) {
              idx_tmp = validIdx.size(0) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r.set_size(nx);
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  r[nx] = b_i + 1;
                  nx++;
                }
              }

              idx_tmp = validIdx.size(0) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r1.set_size(nx);
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  r1[nx] = b_i + 1;
                  nx++;
                }
              }

              idx_tmp = validIdx.size(0) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r2.set_size(nx);
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  r2[nx] = b_i + 1;
                  nx++;
                }
              }

              if (r.size(0) == r2.size(0)) {
                b_idx = static_cast<int>(idx[0]);
                c_idx = static_cast<int>(idx[2]);
                d_idx = static_cast<int>(idx[1]);
                c_this.set_size(r.size(0), 1, BoardCoords.size(2));
                loop_ub = BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i1,k)

                for (int i = 0; i < loop_ub; i++) {
                  k = r.size(0);
                  for (i1 = 0; i1 < k; i1++) {
                    c_this[i1 + c_this.size(0) * i] = (BoardCoords[((r[i1] +
                      BoardCoords.size(0) * (b_idx - 1)) + BoardCoords.size(0) *
                      BoardCoords.size(1) * i) - 1] + BoardCoords[((r1[i1] +
                      BoardCoords.size(0) * (c_idx - 1)) + BoardCoords.size(0) *
                      BoardCoords.size(1) * i) - 1]) - 2.0 * BoardCoords[((r2[i1]
                      + BoardCoords.size(0) * (d_idx - 1)) + BoardCoords.size(0)
                      * BoardCoords.size(1) * i) - 1];
                  }
                }

                b_squeeze(c_this, num);
              } else {
                binary_expand_op(num, this, r, idx, r1, r2);
              }

              idx_tmp = validIdx.size(0) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r3.set_size(nx);
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  r3[nx] = b_i + 1;
                  nx++;
                }
              }

              idx_tmp = validIdx.size(0) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r4.set_size(nx);
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  r4[nx] = b_i + 1;
                  nx++;
                }
              }

              b_idx = static_cast<int>(idx[0]);
              c_idx = static_cast<int>(idx[2]);
              c_this.set_size(r3.size(0), 1, BoardCoords.size(2));
              loop_ub = BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i1,k)

              for (int i = 0; i < loop_ub; i++) {
                k = r3.size(0);
                for (i1 = 0; i1 < k; i1++) {
                  c_this[i1 + c_this.size(0) * i] = BoardCoords[((r3[i1] +
                    BoardCoords.size(0) * (b_idx - 1)) + BoardCoords.size(0) *
                    BoardCoords.size(1) * i) - 1] - BoardCoords[((r4[i1] +
                    BoardCoords.size(0) * (c_idx - 1)) + BoardCoords.size(0) *
                    BoardCoords.size(1) * i) - 1];
                }
              }

              b_squeeze(c_this, b_denom);
              if (num.size(1) > 1) {
                b_r.set_size(num.size(0));
                nx = num.size(0);
                if (static_cast<int>(num.size(0) < 3200)) {
                  for (k = 0; k < nx; k++) {
                    b_r[k] = rt_hypotd_snf(num[k], num[k + num.size(0)]);
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (k = 0; k < nx; k++) {
                    b_r[k] = rt_hypotd_snf(num[k], num[k + num.size(0)]);
                  }
                }

                c_r.set_size(b_denom.size(0));
                nx = b_denom.size(0);
                if (static_cast<int>(b_denom.size(0) < 3200)) {
                  for (k = 0; k < nx; k++) {
                    c_r[k] = rt_hypotd_snf(b_denom[k], b_denom[k + b_denom.size
                      (0)]);
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (k = 0; k < nx; k++) {
                    c_r[k] = rt_hypotd_snf(b_denom[k], b_denom[k + b_denom.size
                      (0)]);
                  }
                }

                if (b_r.size(0) == c_r.size(0)) {
                  loop_ub = b_r.size(0);
                  if (static_cast<int>(b_r.size(0) < 3200)) {
                    for (int i{0}; i < loop_ub; i++) {
                      b_r[i] = b_r[i] / c_r[i];
                    }
                  } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                    for (int i = 0; i < loop_ub; i++) {
                      b_r[i] = b_r[i] / c_r[i];
                    }
                  }

                  newEnergy = std::fmax(oldEnergy, static_cast<float>(::coder::
                    internal::maximum(b_r)));
                } else {
                  newEnergy = binary_expand_op(oldEnergy, b_r, c_r);
                }
              } else {
                newEnergy = std::fmax(oldEnergy, static_cast<float>
                                      (rt_hypotd_snf(num[0], num[1]) /
                  rt_hypotd_snf(b_denom[0], b_denom[1])));
              }
            }

            b_idx = static_cast<int>(idx[0]);
            b_this.set_size(1, BoardIdx.size(0));
            loop_ub = BoardIdx.size(0);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                b_this[i] = (BoardIdx[i + BoardIdx.size(0) * (b_idx - 1)] > 0.0);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                b_this[i] = (BoardIdx[i + BoardIdx.size(0) * (b_idx - 1)] > 0.0);
              }
            }

            Checkerboard::arrayFind(b_this, validNewColIdx);
            if (validNewColIdx.size(1) != 0) {
              int b_loop_ub;
              int e_idx;
              int f_idx;
              idx_tmp = validNewColIdx.size(1);
              c_idx = static_cast<int>(idx[0]);
              d_idx = static_cast<int>(idx[0]);
              b_idx_tmp = static_cast<int>(idx[0]);
              loop_ub = BoardCoords.size(2);
              e_idx = static_cast<int>(idx[0]);
              f_idx = static_cast<int>(idx[0]);
              b_loop_ub = BoardCoords.size(2);
              for (int b_i{0}; b_i < idx_tmp; b_i++) {
                double d;
                d = validNewColIdx[b_i];
                b_num.set_size(1, 1, BoardCoords.size(2));
                for (nx = 0; nx < loop_ub; nx++) {
                  b_num[nx] = (BoardCoords[((static_cast<int>(d) +
                    BoardCoords.size(0) * (c_idx - 1)) + BoardCoords.size(0) *
                    BoardCoords.size(1) * nx) - 1] + BoardCoords[((static_cast<
                    int>(d + 2.0) + BoardCoords.size(0) * (d_idx - 1)) +
                    BoardCoords.size(0) * BoardCoords.size(1) * nx) - 1]) - 2.0 *
                    BoardCoords[((static_cast<int>(d + 1.0) + BoardCoords.size(0)
                                  * (b_idx_tmp - 1)) + BoardCoords.size(0) *
                                 BoardCoords.size(1) * nx) - 1];
                }

                d = validNewColIdx[b_i];
                denom.set_size(1, 1, BoardCoords.size(2));
                for (nx = 0; nx < b_loop_ub; nx++) {
                  denom[nx] = BoardCoords[((static_cast<int>(d) +
                    BoardCoords.size(0) * (e_idx - 1)) + BoardCoords.size(0) *
                    BoardCoords.size(1) * nx) - 1] - BoardCoords[((static_cast<
                    int>(d + 2.0) + BoardCoords.size(0) * (f_idx - 1)) +
                    BoardCoords.size(0) * BoardCoords.size(1) * nx) - 1];
                }

                if (newEnergy != 0.0F) {
                  b_idx = b_num.size(2);
                  nx = denom.size(2);
                  b_r = b_num.reshape(b_idx);
                  c_r = denom.reshape(nx);
                  newEnergy = std::fmax(newEnergy, static_cast<float>(b_norm(b_r)
                    / b_norm(c_r)));
                } else {
                  b_idx = b_num.size(2);
                  nx = denom.size(2);
                  b_r = b_num.reshape(b_idx);
                  c_r = denom.reshape(nx);
                  newEnergy = std::fmax(oldEnergy, static_cast<float>(b_norm(b_r)
                    / b_norm(c_r)));
                }
              }
            }

            if (newEnergy != 0.0F) {
              newEnergy = newEnergy * static_cast<float>(BoardIdx.size(0) *
                BoardIdx.size(1)) - static_cast<float>(BoardIdx.size(0) *
                BoardIdx.size(1));
            } else {
              newEnergy = rtInfF;
            }

            return newEnergy;
          }

          //
          // Arguments    : float oldEnergy
          // Return Type  : float
          //
          float Checkerboard::computeNewEnergyVertical(float oldEnergy) const
          {
            array<double, 3U> b;
            array<double, 3U> denom;
            array<double, 3U> num;
            array<double, 3U> r1;
            array<double, 3U> r3;
            array<double, 3U> r5;
            array<double, 2U> b_denom;
            array<double, 2U> b_num;
            array<double, 2U> validNewRowIdx;
            array<double, 1U> b_r;
            array<double, 1U> c_r;
            array<int, 2U> r;
            array<int, 2U> r2;
            array<int, 2U> r4;
            array<int, 2U> r6;
            array<int, 2U> r7;
            array<boolean_T, 2U> validIdx;
            float newEnergy;
            int end;
            int i2;
            int k;
            int loop_ub;
            int nx;
            boolean_T exitg1;
            boolean_T y;
            validIdx.set_size(1, BoardIdx.size(1));
            loop_ub = BoardIdx.size(1);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                validIdx[i] = ((BoardIdx[BoardIdx.size(0) * i] > 0.0) &&
                               (BoardIdx[BoardIdx.size(0) * i + 1] > 0.0) &&
                               (BoardIdx[BoardIdx.size(0) * i + 2] > 0.0));
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                validIdx[i] = ((BoardIdx[BoardIdx.size(0) * i] > 0.0) &&
                               (BoardIdx[BoardIdx.size(0) * i + 1] > 0.0) &&
                               (BoardIdx[BoardIdx.size(0) * i + 2] > 0.0));
              }
            }

            newEnergy = 0.0F;
            y = false;
            nx = 1;
            exitg1 = false;
            while ((!exitg1) && (nx <= validIdx.size(1))) {
              if (validIdx[nx - 1]) {
                y = true;
                exitg1 = true;
              } else {
                nx++;
              }
            }

            if (y) {
              end = validIdx.size(1) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r.set_size(1, nx);
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  r[nx] = b_i + 1;
                  nx++;
                }
              }

              r1.set_size(1, r.size(1), BoardCoords.size(2));
              loop_ub = BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i2,k)

              for (int i = 0; i < loop_ub; i++) {
                k = r.size(1);
                for (i2 = 0; i2 < k; i2++) {
                  r1[i2 + r1.size(1) * i] = BoardCoords[BoardCoords.size(0) *
                    (r[i2] - 1) + BoardCoords.size(0) * BoardCoords.size(1) * i];
                }
              }

              end = validIdx.size(1) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r2.set_size(1, nx);
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  r2[nx] = b_i + 1;
                  nx++;
                }
              }

              r3.set_size(1, r2.size(1), BoardCoords.size(2));
              loop_ub = BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i2,k)

              for (int i = 0; i < loop_ub; i++) {
                k = r2.size(1);
                for (i2 = 0; i2 < k; i2++) {
                  r3[i2 + r3.size(1) * i] = BoardCoords[(BoardCoords.size(0) *
                    (r2[i2] - 1) + BoardCoords.size(0) * BoardCoords.size(1) * i)
                    + 2];
                }
              }

              end = validIdx.size(1) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r4.set_size(1, nx);
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  r4[nx] = b_i + 1;
                  nx++;
                }
              }

              b.set_size(1, r4.size(1), BoardCoords.size(2));
              loop_ub = BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i2,k)

              for (int i = 0; i < loop_ub; i++) {
                k = r4.size(1);
                for (i2 = 0; i2 < k; i2++) {
                  b[i2 + b.size(1) * i] = BoardCoords[(BoardCoords.size(0) *
                    (r4[i2] - 1) + BoardCoords.size(0) * BoardCoords.size(1) * i)
                    + 1];
                }
              }

              if ((r1.size(1) == b.size(1)) && (r1.size(2) == b.size(2))) {
                r5.set_size(1, r1.size(1), r1.size(2));
                loop_ub = r1.size(1) * r1.size(2);
                if (static_cast<int>(loop_ub < 3200)) {
                  for (int i{0}; i < loop_ub; i++) {
                    r5[i] = (r1[i] + r3[i]) - 2.0 * b[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (int i = 0; i < loop_ub; i++) {
                    r5[i] = (r1[i] + r3[i]) - 2.0 * b[i];
                  }
                }

                squeeze(r5, b_num);
              } else {
                binary_expand_op(b_num, r1, r3, b);
              }

              end = validIdx.size(1) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r6.set_size(1, nx);
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  r6[nx] = b_i + 1;
                  nx++;
                }
              }

              r1.set_size(1, r6.size(1), BoardCoords.size(2));
              loop_ub = BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i2,k)

              for (int i = 0; i < loop_ub; i++) {
                k = r6.size(1);
                for (i2 = 0; i2 < k; i2++) {
                  r1[i2 + r1.size(1) * i] = BoardCoords[BoardCoords.size(0) *
                    (r6[i2] - 1) + BoardCoords.size(0) * BoardCoords.size(1) * i];
                }
              }

              end = validIdx.size(1) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r7.set_size(1, nx);
              nx = 0;
              for (int b_i{0}; b_i <= end; b_i++) {
                if (validIdx[b_i]) {
                  r7[nx] = b_i + 1;
                  nx++;
                }
              }

              r3.set_size(1, r7.size(1), BoardCoords.size(2));
              loop_ub = BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i2,k)

              for (int i = 0; i < loop_ub; i++) {
                k = r7.size(1);
                for (i2 = 0; i2 < k; i2++) {
                  r3[i2 + r3.size(1) * i] = BoardCoords[(BoardCoords.size(0) *
                    (r7[i2] - 1) + BoardCoords.size(0) * BoardCoords.size(1) * i)
                    + 2];
                }
              }

              r5.set_size(1, r1.size(1), r1.size(2));
              loop_ub = r1.size(1) * r1.size(2);
              if (static_cast<int>(loop_ub < 3200)) {
                for (int i{0}; i < loop_ub; i++) {
                  r5[i] = r1[i] - r3[i];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < loop_ub; i++) {
                  r5[i] = r1[i] - r3[i];
                }
              }

              squeeze(r5, b_denom);
              if (b_num.size(1) > 1) {
                b_r.set_size(b_num.size(0));
                nx = b_num.size(0);
                if (static_cast<int>(b_num.size(0) < 3200)) {
                  for (k = 0; k < nx; k++) {
                    b_r[k] = rt_hypotd_snf(b_num[k], b_num[k + b_num.size(0)]);
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (k = 0; k < nx; k++) {
                    b_r[k] = rt_hypotd_snf(b_num[k], b_num[k + b_num.size(0)]);
                  }
                }

                c_r.set_size(b_denom.size(0));
                nx = b_denom.size(0);
                if (static_cast<int>(b_denom.size(0) < 3200)) {
                  for (k = 0; k < nx; k++) {
                    c_r[k] = rt_hypotd_snf(b_denom[k], b_denom[k + b_denom.size
                      (0)]);
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (k = 0; k < nx; k++) {
                    c_r[k] = rt_hypotd_snf(b_denom[k], b_denom[k + b_denom.size
                      (0)]);
                  }
                }

                if (b_r.size(0) == c_r.size(0)) {
                  loop_ub = b_r.size(0);
                  if (static_cast<int>(b_r.size(0) < 3200)) {
                    for (int i{0}; i < loop_ub; i++) {
                      b_r[i] = b_r[i] / c_r[i];
                    }
                  } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                    for (int i = 0; i < loop_ub; i++) {
                      b_r[i] = b_r[i] / c_r[i];
                    }
                  }

                  newEnergy = std::fmax(oldEnergy, static_cast<float>(::coder::
                    internal::maximum(b_r)));
                } else {
                  newEnergy = binary_expand_op(oldEnergy, b_r, c_r);
                }
              } else {
                newEnergy = std::fmax(oldEnergy, static_cast<float>
                                      (rt_hypotd_snf(b_num[0], b_num[1]) /
                  rt_hypotd_snf(b_denom[0], b_denom[1])));
              }
            }

            validIdx.set_size(1, BoardIdx.size(1));
            loop_ub = BoardIdx.size(1);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                validIdx[i] = (BoardIdx[BoardIdx.size(0) * i] > 0.0);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                validIdx[i] = (BoardIdx[BoardIdx.size(0) * i] > 0.0);
              }
            }

            Checkerboard::arrayFind(validIdx, validNewRowIdx);
            if (validNewRowIdx.size(1) != 0) {
              int b_loop_ub;
              int i1;
              i1 = validNewRowIdx.size(1);
              loop_ub = BoardCoords.size(2);
              b_loop_ub = BoardCoords.size(2);
              for (int b_i{0}; b_i < i1; b_i++) {
                double d;
                d = validNewRowIdx[b_i];
                num.set_size(1, 1, BoardCoords.size(2));
                for (nx = 0; nx < loop_ub; nx++) {
                  num[nx] = (BoardCoords[BoardCoords.size(0) * (static_cast<int>
                              (d) - 1) + BoardCoords.size(0) * BoardCoords.size
                             (1) * nx] + BoardCoords[BoardCoords.size(0) * (
                              static_cast<int>(d + 2.0) - 1) + BoardCoords.size
                             (0) * BoardCoords.size(1) * nx]) - 2.0 *
                    BoardCoords[BoardCoords.size(0) * (static_cast<int>(d + 1.0)
                    - 1) + BoardCoords.size(0) * BoardCoords.size(1) * nx];
                }

                d = validNewRowIdx[b_i];
                denom.set_size(1, 1, BoardCoords.size(2));
                for (nx = 0; nx < b_loop_ub; nx++) {
                  denom[nx] = BoardCoords[BoardCoords.size(0) * (static_cast<int>
                    (d) - 1) + BoardCoords.size(0) * BoardCoords.size(1) * nx] -
                    BoardCoords[BoardCoords.size(0) * (static_cast<int>(d + 2.0)
                    - 1) + BoardCoords.size(0) * BoardCoords.size(1) * nx];
                }

                if (newEnergy != 0.0F) {
                  nx = num.size(2);
                  end = denom.size(2);
                  b_r = num.reshape(nx);
                  c_r = denom.reshape(end);
                  newEnergy = std::fmax(newEnergy, static_cast<float>(b_norm(b_r)
                    / b_norm(c_r)));
                } else {
                  nx = num.size(2);
                  end = denom.size(2);
                  b_r = num.reshape(nx);
                  c_r = denom.reshape(end);
                  newEnergy = std::fmax(oldEnergy, static_cast<float>(b_norm(b_r)
                    / b_norm(c_r)));
                }
              }
            }

            if (newEnergy != 0.0F) {
              newEnergy = newEnergy * static_cast<float>(BoardIdx.size(0) *
                BoardIdx.size(1)) - static_cast<float>(BoardIdx.size(0) *
                BoardIdx.size(1));
            } else {
              newEnergy = rtInfF;
            }

            return newEnergy;
          }

          //
          // Arguments    : const ::coder::array<double, 2U> &idx
          //                float oldEnergy
          // Return Type  : float
          //
          float Checkerboard::computeNewEnergyVertical(const ::coder::array<
            double, 2U> &idx, float oldEnergy) const
          {
            array<double, 3U> b;
            array<double, 3U> denom;
            array<double, 3U> num;
            array<double, 3U> r1;
            array<double, 3U> r3;
            array<double, 3U> r5;
            array<double, 2U> b_denom;
            array<double, 2U> b_num;
            array<double, 2U> validNewRowIdx;
            array<double, 1U> b_r;
            array<double, 1U> c_r;
            array<int, 2U> r;
            array<int, 2U> r2;
            array<int, 2U> r4;
            array<int, 2U> r6;
            array<int, 2U> r7;
            array<boolean_T, 2U> validIdx;
            float newEnergy;
            int b_idx;
            int b_idx_tmp;
            int i1;
            int idx_tmp;
            int k;
            int loop_ub;
            int nx;
            boolean_T exitg1;
            boolean_T y;
            nx = static_cast<int>(idx[0]);
            idx_tmp = static_cast<int>(idx[1]);
            b_idx_tmp = static_cast<int>(idx[2]);
            validIdx.set_size(1, BoardIdx.size(1));
            loop_ub = BoardIdx.size(1);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                validIdx[i] = ((BoardIdx[(nx + BoardIdx.size(0) * i) - 1] > 0.0)
                               && (BoardIdx[(idx_tmp + BoardIdx.size(0) * i) - 1]
                                   > 0.0) && (BoardIdx[(b_idx_tmp +
                  BoardIdx.size(0) * i) - 1] > 0.0));
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                validIdx[i] = ((BoardIdx[(nx + BoardIdx.size(0) * i) - 1] > 0.0)
                               && (BoardIdx[(idx_tmp + BoardIdx.size(0) * i) - 1]
                                   > 0.0) && (BoardIdx[(b_idx_tmp +
                  BoardIdx.size(0) * i) - 1] > 0.0));
              }
            }

            newEnergy = 0.0F;
            y = false;
            nx = 1;
            exitg1 = false;
            while ((!exitg1) && (nx <= validIdx.size(1))) {
              if (validIdx[nx - 1]) {
                y = true;
                exitg1 = true;
              } else {
                nx++;
              }
            }

            if (y) {
              idx_tmp = validIdx.size(1) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r.set_size(1, nx);
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  r[nx] = b_i + 1;
                  nx++;
                }
              }

              b_idx = static_cast<int>(idx[0]);
              r1.set_size(1, r.size(1), BoardCoords.size(2));
              loop_ub = BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i1,k)

              for (int i = 0; i < loop_ub; i++) {
                k = r.size(1);
                for (i1 = 0; i1 < k; i1++) {
                  r1[i1 + r1.size(1) * i] = BoardCoords[((b_idx +
                    BoardCoords.size(0) * (r[i1] - 1)) + BoardCoords.size(0) *
                    BoardCoords.size(1) * i) - 1];
                }
              }

              idx_tmp = validIdx.size(1) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r2.set_size(1, nx);
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  r2[nx] = b_i + 1;
                  nx++;
                }
              }

              b_idx = static_cast<int>(idx[2]);
              r3.set_size(1, r2.size(1), BoardCoords.size(2));
              loop_ub = BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i1,k)

              for (int i = 0; i < loop_ub; i++) {
                k = r2.size(1);
                for (i1 = 0; i1 < k; i1++) {
                  r3[i1 + r3.size(1) * i] = BoardCoords[((b_idx +
                    BoardCoords.size(0) * (r2[i1] - 1)) + BoardCoords.size(0) *
                    BoardCoords.size(1) * i) - 1];
                }
              }

              idx_tmp = validIdx.size(1) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r4.set_size(1, nx);
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  r4[nx] = b_i + 1;
                  nx++;
                }
              }

              b_idx = static_cast<int>(idx[1]);
              b.set_size(1, r4.size(1), BoardCoords.size(2));
              loop_ub = BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i1,k)

              for (int i = 0; i < loop_ub; i++) {
                k = r4.size(1);
                for (i1 = 0; i1 < k; i1++) {
                  b[i1 + b.size(1) * i] = BoardCoords[((b_idx + BoardCoords.size
                    (0) * (r4[i1] - 1)) + BoardCoords.size(0) * BoardCoords.size
                    (1) * i) - 1];
                }
              }

              if ((r1.size(1) == b.size(1)) && (r1.size(2) == b.size(2))) {
                r5.set_size(1, r1.size(1), r1.size(2));
                loop_ub = r1.size(1) * r1.size(2);
                if (static_cast<int>(loop_ub < 3200)) {
                  for (int i{0}; i < loop_ub; i++) {
                    r5[i] = (r1[i] + r3[i]) - 2.0 * b[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (int i = 0; i < loop_ub; i++) {
                    r5[i] = (r1[i] + r3[i]) - 2.0 * b[i];
                  }
                }

                squeeze(r5, b_num);
              } else {
                binary_expand_op(b_num, r1, r3, b);
              }

              idx_tmp = validIdx.size(1) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r6.set_size(1, nx);
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  r6[nx] = b_i + 1;
                  nx++;
                }
              }

              b_idx = static_cast<int>(idx[0]);
              r1.set_size(1, r6.size(1), BoardCoords.size(2));
              loop_ub = BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i1,k)

              for (int i = 0; i < loop_ub; i++) {
                k = r6.size(1);
                for (i1 = 0; i1 < k; i1++) {
                  r1[i1 + r1.size(1) * i] = BoardCoords[((b_idx +
                    BoardCoords.size(0) * (r6[i1] - 1)) + BoardCoords.size(0) *
                    BoardCoords.size(1) * i) - 1];
                }
              }

              idx_tmp = validIdx.size(1) - 1;
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  nx++;
                }
              }

              r7.set_size(1, nx);
              nx = 0;
              for (int b_i{0}; b_i <= idx_tmp; b_i++) {
                if (validIdx[b_i]) {
                  r7[nx] = b_i + 1;
                  nx++;
                }
              }

              b_idx = static_cast<int>(idx[2]);
              r3.set_size(1, r7.size(1), BoardCoords.size(2));
              loop_ub = BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i1,k)

              for (int i = 0; i < loop_ub; i++) {
                k = r7.size(1);
                for (i1 = 0; i1 < k; i1++) {
                  r3[i1 + r3.size(1) * i] = BoardCoords[((b_idx +
                    BoardCoords.size(0) * (r7[i1] - 1)) + BoardCoords.size(0) *
                    BoardCoords.size(1) * i) - 1];
                }
              }

              r5.set_size(1, r1.size(1), r1.size(2));
              loop_ub = r1.size(1) * r1.size(2);
              if (static_cast<int>(loop_ub < 3200)) {
                for (int i{0}; i < loop_ub; i++) {
                  r5[i] = r1[i] - r3[i];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < loop_ub; i++) {
                  r5[i] = r1[i] - r3[i];
                }
              }

              squeeze(r5, b_denom);
              if (b_num.size(1) > 1) {
                b_r.set_size(b_num.size(0));
                nx = b_num.size(0);
                if (static_cast<int>(b_num.size(0) < 3200)) {
                  for (k = 0; k < nx; k++) {
                    b_r[k] = rt_hypotd_snf(b_num[k], b_num[k + b_num.size(0)]);
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (k = 0; k < nx; k++) {
                    b_r[k] = rt_hypotd_snf(b_num[k], b_num[k + b_num.size(0)]);
                  }
                }

                c_r.set_size(b_denom.size(0));
                nx = b_denom.size(0);
                if (static_cast<int>(b_denom.size(0) < 3200)) {
                  for (k = 0; k < nx; k++) {
                    c_r[k] = rt_hypotd_snf(b_denom[k], b_denom[k + b_denom.size
                      (0)]);
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (k = 0; k < nx; k++) {
                    c_r[k] = rt_hypotd_snf(b_denom[k], b_denom[k + b_denom.size
                      (0)]);
                  }
                }

                if (b_r.size(0) == c_r.size(0)) {
                  loop_ub = b_r.size(0);
                  if (static_cast<int>(b_r.size(0) < 3200)) {
                    for (int i{0}; i < loop_ub; i++) {
                      b_r[i] = b_r[i] / c_r[i];
                    }
                  } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                    for (int i = 0; i < loop_ub; i++) {
                      b_r[i] = b_r[i] / c_r[i];
                    }
                  }

                  newEnergy = std::fmax(oldEnergy, static_cast<float>(::coder::
                    internal::maximum(b_r)));
                } else {
                  newEnergy = binary_expand_op(oldEnergy, b_r, c_r);
                }
              } else {
                newEnergy = std::fmax(oldEnergy, static_cast<float>
                                      (rt_hypotd_snf(b_num[0], b_num[1]) /
                  rt_hypotd_snf(b_denom[0], b_denom[1])));
              }
            }

            b_idx = static_cast<int>(idx[0]);
            validIdx.set_size(1, BoardIdx.size(1));
            loop_ub = BoardIdx.size(1);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                validIdx[i] = (BoardIdx[(b_idx + BoardIdx.size(0) * i) - 1] >
                               0.0);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                validIdx[i] = (BoardIdx[(b_idx + BoardIdx.size(0) * i) - 1] >
                               0.0);
              }
            }

            Checkerboard::arrayFind(validIdx, validNewRowIdx);
            if (validNewRowIdx.size(1) != 0) {
              int b_loop_ub;
              int c_idx;
              int d_idx;
              int e_idx;
              int f_idx;
              idx_tmp = validNewRowIdx.size(1);
              b_idx_tmp = static_cast<int>(idx[0]);
              c_idx = static_cast<int>(idx[0]);
              d_idx = static_cast<int>(idx[0]);
              loop_ub = BoardCoords.size(2);
              e_idx = static_cast<int>(idx[0]);
              f_idx = static_cast<int>(idx[0]);
              b_loop_ub = BoardCoords.size(2);
              for (int b_i{0}; b_i < idx_tmp; b_i++) {
                double d;
                d = validNewRowIdx[b_i];
                num.set_size(1, 1, BoardCoords.size(2));
                for (nx = 0; nx < loop_ub; nx++) {
                  num[nx] = (BoardCoords[((b_idx_tmp + BoardCoords.size(0) * (
                    static_cast<int>(d) - 1)) + BoardCoords.size(0) *
                              BoardCoords.size(1) * nx) - 1] + BoardCoords
                             [((c_idx + BoardCoords.size(0) * (static_cast<int>
                    (d + 2.0) - 1)) + BoardCoords.size(0) * BoardCoords.size(1) *
                               nx) - 1]) - 2.0 * BoardCoords[((d_idx +
                    BoardCoords.size(0) * (static_cast<int>(d + 1.0) - 1)) +
                    BoardCoords.size(0) * BoardCoords.size(1) * nx) - 1];
                }

                d = validNewRowIdx[b_i];
                denom.set_size(1, 1, BoardCoords.size(2));
                for (nx = 0; nx < b_loop_ub; nx++) {
                  denom[nx] = BoardCoords[((e_idx + BoardCoords.size(0) * (
                    static_cast<int>(d) - 1)) + BoardCoords.size(0) *
                    BoardCoords.size(1) * nx) - 1] - BoardCoords[((f_idx +
                    BoardCoords.size(0) * (static_cast<int>(d + 2.0) - 1)) +
                    BoardCoords.size(0) * BoardCoords.size(1) * nx) - 1];
                }

                if (newEnergy != 0.0F) {
                  b_idx = num.size(2);
                  nx = denom.size(2);
                  b_r = num.reshape(b_idx);
                  c_r = denom.reshape(nx);
                  newEnergy = std::fmax(newEnergy, static_cast<float>(b_norm(b_r)
                    / b_norm(c_r)));
                } else {
                  b_idx = num.size(2);
                  nx = denom.size(2);
                  b_r = num.reshape(b_idx);
                  c_r = denom.reshape(nx);
                  newEnergy = std::fmax(oldEnergy, static_cast<float>(b_norm(b_r)
                    / b_norm(c_r)));
                }
              }
            }

            if (newEnergy != 0.0F) {
              newEnergy = newEnergy * static_cast<float>(BoardIdx.size(0) *
                BoardIdx.size(1)) - static_cast<float>(BoardIdx.size(0) *
                BoardIdx.size(1));
            } else {
              newEnergy = rtInfF;
            }

            return newEnergy;
          }

          //
          // Arguments    : ::coder::array<double, 2U> &newIndices
          // Return Type  : void
          //
          void Checkerboard::d_fitPolynomialIndices(::coder::array<double, 2U>
            &newIndices) const
          {
            array<double, 2U> b_index;
            array<double, 2U> b_this;
            array<double, 2U> removedIdx;
            array<int, 2U> validIdx;
            array<int, 1U> currCurve_tmp;
            double currCurve_data[5];
            double coordsToUse[2];
            int b_i;
            int loop_ub;
            int n;
            b_findIndependentVar(coordsToUse);
            newIndices.set_size(1, BoardCoords.size(0));
            loop_ub = BoardCoords.size(0);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                newIndices[i] = 0.0;
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                newIndices[i] = 0.0;
              }
            }

            removedIdx.set_size(1, 0);
            b_i = BoardCoords.size(0);
            for (int j{0}; j < b_i; j++) {
              int i1;
              b_index.set_size(1, BoardCoords.size(1));
              loop_ub = BoardCoords.size(1);
              for (i1 = 0; i1 < loop_ub; i1++) {
                b_index[i1] = BoardCoords[(j + BoardCoords.size(0) * i1) +
                  BoardCoords.size(0) * BoardCoords.size(1) * (static_cast<int>
                  (coordsToUse[0]) - 1)];
              }

              eml_find(b_index, validIdx);
              if (validIdx.size(1) >= 2) {
                double coordDist;
                double currCoord;
                double currRad;
                double refCoordValue;
                int currCurve_size[2];
                int coordDist_tmp;
                int i2;
                boolean_T exitg1;
                coordDist_tmp = validIdx[0];
                coordDist = (BoardCoords[(j + BoardCoords.size(0) * (validIdx[0]
                  - 1)) + BoardCoords.size(0) * BoardCoords.size(1) * (
                  static_cast<int>(coordsToUse[0]) - 1)] - BoardCoords[(j +
                  BoardCoords.size(0) * (validIdx[1] - 1)) + BoardCoords.size(0)
                             * BoardCoords.size(1) * (static_cast<int>
                  (coordsToUse[0]) - 1)]) / (static_cast<double>(validIdx[1]) -
                  static_cast<double>(coordDist_tmp));
                n = 0;
                i1 = validIdx.size(1);
                currCurve_tmp.set_size(validIdx.size(1));
                for (loop_ub = 0; loop_ub < i1; loop_ub++) {
                  i2 = validIdx[loop_ub];
                  if (i2 != 0) {
                    n++;
                  }

                  currCurve_tmp[loop_ub] = i2;
                }

                b_index.set_size(1, currCurve_tmp.size(0));
                loop_ub = currCurve_tmp.size(0);
                for (i1 = 0; i1 < loop_ub; i1++) {
                  b_index[i1] = BoardCoords[(j + BoardCoords.size(0) *
                    (currCurve_tmp[i1] - 1)) + BoardCoords.size(0) *
                    BoardCoords.size(1) * (static_cast<int>(coordsToUse[0]) - 1)];
                }

                b_this.set_size(1, currCurve_tmp.size(0));
                loop_ub = currCurve_tmp.size(0);
                for (i1 = 0; i1 < loop_ub; i1++) {
                  b_this[i1] = BoardCoords[(j + BoardCoords.size(0) *
                    (currCurve_tmp[i1] - 1)) + BoardCoords.size(0) *
                    BoardCoords.size(1) * (static_cast<int>(coordsToUse[1]) - 1)];
                }

                if (n > 5) {
                  i1 = 4;
                } else {
                  i1 = 2;
                }

                polyfit(b_index, b_this, static_cast<double>(i1), currCurve_data,
                        currCurve_size);
                currRad = coordDist / 4.0;
                refCoordValue = BoardCoords[(j + BoardCoords.size(0) *
                  (validIdx[0] - 1)) + BoardCoords.size(0) * BoardCoords.size(1)
                  * (static_cast<int>(coordsToUse[0]) - 1)];
                currCoord = currRad + refCoordValue;
                exitg1 = false;
                while ((!exitg1) && (std::abs(currCoord - refCoordValue) <
                                     static_cast<double>(coordDist_tmp) * 1.5 *
                                     std::abs(coordDist))) {
                  double currPt[2];
                  boolean_T exitg2;
                  boolean_T p;
                  p = true;
                  loop_ub = 0;
                  exitg2 = false;
                  while ((!exitg2) && (loop_ub < 2)) {
                    if (!(coordsToUse[loop_ub] == static_cast<double>(loop_ub) +
                          1.0)) {
                      p = false;
                      exitg2 = true;
                    } else {
                      loop_ub++;
                    }
                  }

                  if (p) {
                    double y;
                    y = currCurve_data[0];
                    i1 = currCurve_size[1];
                    for (loop_ub = 0; loop_ub <= i1 - 2; loop_ub++) {
                      y = currCoord * y + currCurve_data[loop_ub + 1];
                    }

                    currPt[0] = currCoord;
                    currPt[1] = y;
                  } else {
                    double y;
                    y = currCurve_data[0];
                    i1 = currCurve_size[1];
                    for (loop_ub = 0; loop_ub <= i1 - 2; loop_ub++) {
                      y = currCoord * y + currCurve_data[loop_ub + 1];
                    }

                    currPt[0] = y;
                    currPt[1] = currCoord;
                  }

                  findClosestOnCurve(currPt, std::abs(currRad), currCurve_data,
                                     currCurve_size, coordsToUse, removedIdx,
                                     b_index);
                  if (b_index.size(1) != 0) {
                    newIndices[j] = b_index[0];
                    i1 = removedIdx.size(1);
                    loop_ub = b_index.size(1);
                    removedIdx.set_size(removedIdx.size(0), removedIdx.size(1) +
                                        b_index.size(1));
                    for (i2 = 0; i2 < loop_ub; i2++) {
                      removedIdx[i1 + i2] = b_index[i2];
                    }

                    exitg1 = true;
                  } else {
                    currCoord += currRad;
                  }
                }
              }
            }

            n = 0;
            b_i = newIndices.size(1);
            for (loop_ub = 0; loop_ub < b_i; loop_ub++) {
              if (newIndices[loop_ub] != 0.0) {
                n++;
              }
            }

            if (n < 4) {
              loop_ub = newIndices.size(1);
              if (static_cast<int>(newIndices.size(1) < 3200)) {
                for (int i{0}; i < loop_ub; i++) {
                  if (newIndices[i] > 0.0) {
                    newIndices[i] = 0.0;
                  }
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < loop_ub; i++) {
                  if (newIndices[i] > 0.0) {
                    newIndices[i] = 0.0;
                  }
                }
              }
            }
          }

          //
          // Arguments    : const ::coder::array<double, 2U> &indices
          //                ::coder::array<double, 2U> &newBoard
          //                ::coder::array<double, 3U> &newBoardCoords
          // Return Type  : void
          //
          void Checkerboard::expandBoardDown(const ::coder::array<double, 2U>
            &indices, ::coder::array<double, 2U> &newBoard, ::coder::array<
            double, 3U> &newBoardCoords) const
          {
            array<double, 2U> r1;
            array<int, 2U> r;
            array<int, 2U> r2;
            int b_loop_ub;
            int b_this;
            int c_loop_ub;
            int end;
            int i1;
            int i4;
            int loop_ub;
            newBoard.set_size(BoardIdx.size(0) + 1, BoardIdx.size(1));
            loop_ub = (BoardIdx.size(0) + 1) * BoardIdx.size(1);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                newBoard[i] = 0.0;
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                newBoard[i] = 0.0;
              }
            }

            b_this = BoardIdx.size(0);
            loop_ub = indices.size(1);
            if (static_cast<int>(indices.size(1) < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                newBoard[b_this + newBoard.size(0) * i] = indices[i];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                newBoard[b_this + newBoard.size(0) * i] = indices[i];
              }
            }

            loop_ub = BoardIdx.size(1);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i1,b_loop_ub)

            for (int i = 0; i < loop_ub; i++) {
              b_loop_ub = BoardIdx.size(0);
              for (i1 = 0; i1 < b_loop_ub; i1++) {
                newBoard[i1 + newBoard.size(0) * i] = BoardIdx[i1 +
                  BoardIdx.size(0) * i];
              }
            }

            newBoardCoords.set_size(BoardCoords.size(0) + 1, BoardCoords.size(1),
              BoardCoords.size(2));
            loop_ub = (BoardCoords.size(0) + 1) * BoardCoords.size(1) *
              BoardCoords.size(2);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                newBoardCoords[i] = 0.0;
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                newBoardCoords[i] = 0.0;
              }
            }

            end = indices.size(1) - 1;
            loop_ub = 0;
            for (b_this = 0; b_this <= end; b_this++) {
              if (indices[b_this] > 0.0) {
                loop_ub++;
              }
            }

            r.set_size(1, loop_ub);
            loop_ub = 0;
            for (b_this = 0; b_this <= end; b_this++) {
              if (indices[b_this] > 0.0) {
                r[loop_ub] = b_this + 1;
                loop_ub++;
              }
            }

            r1.set_size(r.size(1), Points.size(1));
            loop_ub = Points.size(1);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i1,b_loop_ub)

            for (int i = 0; i < loop_ub; i++) {
              b_loop_ub = r.size(1);
              for (i1 = 0; i1 < b_loop_ub; i1++) {
                r1[i1 + r1.size(0) * i] = Points[(static_cast<int>(indices[r[i1]
                  - 1]) + Points.size(0) * i) - 1];
              }
            }

            end = indices.size(1) - 1;
            loop_ub = 0;
            for (b_this = 0; b_this <= end; b_this++) {
              if (indices[b_this] > 0.0) {
                loop_ub++;
              }
            }

            r2.set_size(1, loop_ub);
            loop_ub = 0;
            for (b_this = 0; b_this <= end; b_this++) {
              if (indices[b_this] > 0.0) {
                r2[loop_ub] = b_this + 1;
                loop_ub++;
              }
            }

            loop_ub = r2.size(1);
            end = BoardCoords.size(2);
            b_this = BoardCoords.size(0);
            for (int i2{0}; i2 < end; i2++) {
              for (int i3{0}; i3 < loop_ub; i3++) {
                newBoardCoords[(b_this + newBoardCoords.size(0) * (r2[i3] - 1))
                  + newBoardCoords.size(0) * newBoardCoords.size(1) * i2] =
                  r1[i3 + loop_ub * i2];
              }
            }

            loop_ub = BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i1,b_loop_ub,i4,c_loop_ub)

            for (int i = 0; i < loop_ub; i++) {
              c_loop_ub = BoardCoords.size(1);
              for (i4 = 0; i4 < c_loop_ub; i4++) {
                b_loop_ub = BoardCoords.size(0);
                for (i1 = 0; i1 < b_loop_ub; i1++) {
                  newBoardCoords[(i1 + newBoardCoords.size(0) * i4) +
                    newBoardCoords.size(0) * newBoardCoords.size(1) * i] =
                    BoardCoords[(i1 + BoardCoords.size(0) * i4) +
                    BoardCoords.size(0) * BoardCoords.size(1) * i];
                }
              }
            }
          }

          //
          // Arguments    : const ::coder::array<double, 2U> &indices
          //                ::coder::array<double, 2U> &newBoard
          //                ::coder::array<double, 3U> &newBoardCoords
          // Return Type  : void
          //
          void Checkerboard::expandBoardLeft(const ::coder::array<double, 2U>
            &indices, ::coder::array<double, 2U> &newBoard, ::coder::array<
            double, 3U> &newBoardCoords) const
          {
            array<double, 2U> r1;
            array<int, 2U> r;
            array<int, 2U> r2;
            int b_loop_ub;
            int c_loop_ub;
            int end;
            int i1;
            int i2;
            int i3;
            int loop_ub;
            newBoard.set_size(BoardIdx.size(0), BoardIdx.size(1) + 1);
            loop_ub = BoardIdx.size(0) * (BoardIdx.size(1) + 1);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                newBoard[i] = 0.0;
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                newBoard[i] = 0.0;
              }
            }

            loop_ub = BoardIdx.size(0);
            if (static_cast<int>(BoardIdx.size(0) < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                newBoard[i] = indices[i];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                newBoard[i] = indices[i];
              }
            }

            i1 = (newBoard.size(1) >= 2);
            loop_ub = BoardIdx.size(1);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i2,b_loop_ub)

            for (int i = 0; i < loop_ub; i++) {
              b_loop_ub = BoardIdx.size(0);
              for (i2 = 0; i2 < b_loop_ub; i2++) {
                newBoard[i2 + newBoard.size(0) * (i1 + i)] = BoardIdx[i2 +
                  BoardIdx.size(0) * i];
              }
            }

            newBoardCoords.set_size(BoardCoords.size(0), BoardCoords.size(1) + 1,
              BoardCoords.size(2));
            loop_ub = BoardCoords.size(0) * (BoardCoords.size(1) + 1) *
              BoardCoords.size(2);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                newBoardCoords[i] = 0.0;
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                newBoardCoords[i] = 0.0;
              }
            }

            end = indices.size(1) - 1;
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (indices[b_i] > 0.0) {
                loop_ub++;
              }
            }

            r.set_size(1, loop_ub);
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (indices[b_i] > 0.0) {
                r[loop_ub] = b_i + 1;
                loop_ub++;
              }
            }

            r1.set_size(r.size(1), Points.size(1));
            loop_ub = Points.size(1);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i2,b_loop_ub)

            for (int i = 0; i < loop_ub; i++) {
              b_loop_ub = r.size(1);
              for (i2 = 0; i2 < b_loop_ub; i2++) {
                r1[i2 + r1.size(0) * i] = Points[(static_cast<int>(indices[r[i2]
                  - 1]) + Points.size(0) * i) - 1];
              }
            }

            end = indices.size(1) - 1;
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (indices[b_i] > 0.0) {
                loop_ub++;
              }
            }

            r2.set_size(1, loop_ub);
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (indices[b_i] > 0.0) {
                r2[loop_ub] = b_i + 1;
                loop_ub++;
              }
            }

            loop_ub = r2.size(1);
            end = BoardCoords.size(2);
            for (i1 = 0; i1 < end; i1++) {
              for (int b_i{0}; b_i < loop_ub; b_i++) {
                newBoardCoords[(r2[b_i] + newBoardCoords.size(0) *
                                newBoardCoords.size(1) * i1) - 1] = r1[b_i +
                  loop_ub * i1];
              }
            }

            i1 = (newBoardCoords.size(1) >= 2);
            loop_ub = BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i2,b_loop_ub,i3,c_loop_ub)

            for (int i = 0; i < loop_ub; i++) {
              c_loop_ub = BoardCoords.size(1);
              for (i3 = 0; i3 < c_loop_ub; i3++) {
                b_loop_ub = BoardCoords.size(0);
                for (i2 = 0; i2 < b_loop_ub; i2++) {
                  newBoardCoords[(i2 + newBoardCoords.size(0) * (i1 + i3)) +
                    newBoardCoords.size(0) * newBoardCoords.size(1) * i] =
                    BoardCoords[(i2 + BoardCoords.size(0) * i3) +
                    BoardCoords.size(0) * BoardCoords.size(1) * i];
                }
              }
            }
          }

          //
          // Arguments    : const ::coder::array<double, 2U> &indices
          //                ::coder::array<double, 2U> &newBoard
          //                ::coder::array<double, 3U> &newBoardCoords
          // Return Type  : void
          //
          void Checkerboard::expandBoardRight(const ::coder::array<double, 2U>
            &indices, ::coder::array<double, 2U> &newBoard, ::coder::array<
            double, 3U> &newBoardCoords) const
          {
            array<double, 2U> r1;
            array<int, 2U> r;
            array<int, 2U> r2;
            int b_loop_ub;
            int b_this;
            int c_loop_ub;
            int end;
            int i1;
            int i4;
            int loop_ub;
            newBoard.set_size(BoardIdx.size(0), BoardIdx.size(1) + 1);
            loop_ub = BoardIdx.size(0) * (BoardIdx.size(1) + 1);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                newBoard[i] = 0.0;
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                newBoard[i] = 0.0;
              }
            }

            b_this = BoardIdx.size(1);
            loop_ub = BoardIdx.size(0);
            if (static_cast<int>(BoardIdx.size(0) < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                newBoard[i + newBoard.size(0) * b_this] = indices[i];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                newBoard[i + newBoard.size(0) * b_this] = indices[i];
              }
            }

            loop_ub = BoardIdx.size(1);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i1,b_loop_ub)

            for (int i = 0; i < loop_ub; i++) {
              b_loop_ub = BoardIdx.size(0);
              for (i1 = 0; i1 < b_loop_ub; i1++) {
                newBoard[i1 + newBoard.size(0) * i] = BoardIdx[i1 +
                  BoardIdx.size(0) * i];
              }
            }

            newBoardCoords.set_size(BoardCoords.size(0), BoardCoords.size(1) + 1,
              BoardCoords.size(2));
            loop_ub = BoardCoords.size(0) * (BoardCoords.size(1) + 1) *
              BoardCoords.size(2);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                newBoardCoords[i] = 0.0;
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                newBoardCoords[i] = 0.0;
              }
            }

            end = indices.size(1) - 1;
            loop_ub = 0;
            for (b_this = 0; b_this <= end; b_this++) {
              if (indices[b_this] > 0.0) {
                loop_ub++;
              }
            }

            r.set_size(1, loop_ub);
            loop_ub = 0;
            for (b_this = 0; b_this <= end; b_this++) {
              if (indices[b_this] > 0.0) {
                r[loop_ub] = b_this + 1;
                loop_ub++;
              }
            }

            r1.set_size(r.size(1), Points.size(1));
            loop_ub = Points.size(1);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i1,b_loop_ub)

            for (int i = 0; i < loop_ub; i++) {
              b_loop_ub = r.size(1);
              for (i1 = 0; i1 < b_loop_ub; i1++) {
                r1[i1 + r1.size(0) * i] = Points[(static_cast<int>(indices[r[i1]
                  - 1]) + Points.size(0) * i) - 1];
              }
            }

            end = indices.size(1) - 1;
            loop_ub = 0;
            for (b_this = 0; b_this <= end; b_this++) {
              if (indices[b_this] > 0.0) {
                loop_ub++;
              }
            }

            r2.set_size(1, loop_ub);
            loop_ub = 0;
            for (b_this = 0; b_this <= end; b_this++) {
              if (indices[b_this] > 0.0) {
                r2[loop_ub] = b_this + 1;
                loop_ub++;
              }
            }

            loop_ub = r2.size(1);
            end = BoardCoords.size(2);
            b_this = BoardCoords.size(1);
            for (int i2{0}; i2 < end; i2++) {
              for (int i3{0}; i3 < loop_ub; i3++) {
                newBoardCoords[((r2[i3] + newBoardCoords.size(0) * b_this) +
                                newBoardCoords.size(0) * newBoardCoords.size(1) *
                                i2) - 1] = r1[i3 + loop_ub * i2];
              }
            }

            loop_ub = BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i1,b_loop_ub,i4,c_loop_ub)

            for (int i = 0; i < loop_ub; i++) {
              c_loop_ub = BoardCoords.size(1);
              for (i4 = 0; i4 < c_loop_ub; i4++) {
                b_loop_ub = BoardCoords.size(0);
                for (i1 = 0; i1 < b_loop_ub; i1++) {
                  newBoardCoords[(i1 + newBoardCoords.size(0) * i4) +
                    newBoardCoords.size(0) * newBoardCoords.size(1) * i] =
                    BoardCoords[(i1 + BoardCoords.size(0) * i4) +
                    BoardCoords.size(0) * BoardCoords.size(1) * i];
                }
              }
            }
          }

          //
          // Arguments    : const ::coder::array<double, 2U> &indices
          //                ::coder::array<double, 2U> &newBoard
          //                ::coder::array<double, 3U> &newBoardCoords
          // Return Type  : void
          //
          void Checkerboard::expandBoardUp(const ::coder::array<double, 2U>
            &indices, ::coder::array<double, 2U> &newBoard, ::coder::array<
            double, 3U> &newBoardCoords) const
          {
            array<double, 2U> r1;
            array<int, 2U> r;
            array<int, 2U> r2;
            int b_loop_ub;
            int c_loop_ub;
            int end;
            int i1;
            int i2;
            int i3;
            int loop_ub;
            newBoard.set_size(BoardIdx.size(0) + 1, BoardIdx.size(1));
            loop_ub = (BoardIdx.size(0) + 1) * BoardIdx.size(1);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                newBoard[i] = 0.0;
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                newBoard[i] = 0.0;
              }
            }

            loop_ub = indices.size(1);
            if (static_cast<int>(indices.size(1) < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                newBoard[newBoard.size(0) * i] = indices[i];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                newBoard[newBoard.size(0) * i] = indices[i];
              }
            }

            i1 = (newBoard.size(0) >= 2);
            loop_ub = BoardIdx.size(1);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i2,b_loop_ub)

            for (int i = 0; i < loop_ub; i++) {
              b_loop_ub = BoardIdx.size(0);
              for (i2 = 0; i2 < b_loop_ub; i2++) {
                newBoard[(i1 + i2) + newBoard.size(0) * i] = BoardIdx[i2 +
                  BoardIdx.size(0) * i];
              }
            }

            newBoardCoords.set_size(BoardCoords.size(0) + 1, BoardCoords.size(1),
              BoardCoords.size(2));
            loop_ub = (BoardCoords.size(0) + 1) * BoardCoords.size(1) *
              BoardCoords.size(2);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                newBoardCoords[i] = 0.0;
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                newBoardCoords[i] = 0.0;
              }
            }

            end = indices.size(1) - 1;
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (indices[b_i] > 0.0) {
                loop_ub++;
              }
            }

            r.set_size(1, loop_ub);
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (indices[b_i] > 0.0) {
                r[loop_ub] = b_i + 1;
                loop_ub++;
              }
            }

            r1.set_size(r.size(1), Points.size(1));
            loop_ub = Points.size(1);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i2,b_loop_ub)

            for (int i = 0; i < loop_ub; i++) {
              b_loop_ub = r.size(1);
              for (i2 = 0; i2 < b_loop_ub; i2++) {
                r1[i2 + r1.size(0) * i] = Points[(static_cast<int>(indices[r[i2]
                  - 1]) + Points.size(0) * i) - 1];
              }
            }

            end = indices.size(1) - 1;
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (indices[b_i] > 0.0) {
                loop_ub++;
              }
            }

            r2.set_size(1, loop_ub);
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (indices[b_i] > 0.0) {
                r2[loop_ub] = b_i + 1;
                loop_ub++;
              }
            }

            loop_ub = r2.size(1);
            end = BoardCoords.size(2);
            for (i1 = 0; i1 < end; i1++) {
              for (int b_i{0}; b_i < loop_ub; b_i++) {
                newBoardCoords[newBoardCoords.size(0) * (r2[b_i] - 1) +
                  newBoardCoords.size(0) * newBoardCoords.size(1) * i1] = r1[b_i
                  + loop_ub * i1];
              }
            }

            i1 = (newBoardCoords.size(0) >= 2);
            loop_ub = BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i2,b_loop_ub,i3,c_loop_ub)

            for (int i = 0; i < loop_ub; i++) {
              c_loop_ub = BoardCoords.size(1);
              for (i3 = 0; i3 < c_loop_ub; i3++) {
                b_loop_ub = BoardCoords.size(0);
                for (i2 = 0; i2 < b_loop_ub; i2++) {
                  newBoardCoords[((i1 + i2) + newBoardCoords.size(0) * i3) +
                    newBoardCoords.size(0) * newBoardCoords.size(1) * i] =
                    BoardCoords[(i2 + BoardCoords.size(0) * i3) +
                    BoardCoords.size(0) * BoardCoords.size(1) * i];
                }
              }
            }
          }

          //
          // Arguments    : const ::coder::array<double, 2U> &predictedPoints
          //                ::coder::array<double, 2U> &indices
          // Return Type  : void
          //
          void Checkerboard::findClosestIndices(const ::coder::array<double, 2U>
            &predictedPoints, ::coder::array<double, 2U> &indices) const
          {
            array<double, 2U> remIdx;
            array<double, 2U> y;
            array<double, 1U> b_this;
            array<float, 2U> c_this;
            array<float, 2U> diffs;
            array<float, 1U> dists;
            array<int, 2U> r1;
            array<int, 1U> validPredictions;
            array<boolean_T, 2U> distIdx;
            array<boolean_T, 1U> r;
            float ex;
            int b_loop_ub;
            int c_loop_ub;
            int loop_ub;
            int nx;
            indices.set_size(1, predictedPoints.size(0));
            loop_ub = predictedPoints.size(0);
            if (static_cast<int>(predictedPoints.size(0) < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                indices[i] = 0.0;
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                indices[i] = 0.0;
              }
            }

            nx = Points.size(0);
            if (nx < 1) {
              y.set_size(1, 0);
            } else {
              y.set_size(1, nx);
              loop_ub = nx - 1;
              if (static_cast<int>(nx < 3200)) {
                for (int i{0}; i <= loop_ub; i++) {
                  y[i] = static_cast<double>(i) + 1.0;
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i <= loop_ub; i++) {
                  y[i] = static_cast<double>(i) + 1.0;
                }
              }
            }

            nx = BoardIdx.size(0) * BoardIdx.size(1);
            b_this = BoardIdx.reshape(nx);
            do_vectors(y, b_this, remIdx, validPredictions, &nx);
            if (remIdx.size(1) != 0) {
              int i1;
              r.set_size(predictedPoints.size(0));
              loop_ub = predictedPoints.size(0);
              if (static_cast<int>(predictedPoints.size(0) < 3200)) {
                for (int i{0}; i < loop_ub; i++) {
                  r[i] = std::isnan(predictedPoints[i]);
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < loop_ub; i++) {
                  r[i] = std::isnan(predictedPoints[i]);
                }
              }

              loop_ub = r.size(0);
              if (static_cast<int>(r.size(0) < 3200)) {
                for (int i{0}; i < loop_ub; i++) {
                  r[i] = !r[i];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < loop_ub; i++) {
                  r[i] = !r[i];
                }
              }

              b_eml_find(r, validPredictions);
              i1 = validPredictions.size(0);
              if (validPredictions.size(0) - 1 >= 0) {
                b_loop_ub = Points.size(1);
                c_loop_ub = predictedPoints.size(1);
              }

              for (int b_i{0}; b_i < i1; b_i++) {
                c_this.set_size(remIdx.size(1), Points.size(1));
                for (nx = 0; nx < b_loop_ub; nx++) {
                  loop_ub = remIdx.size(1);
                  for (int c_i{0}; c_i < loop_ub; c_i++) {
                    c_this[c_i + c_this.size(0) * nx] = Points[(static_cast<int>
                      (remIdx[c_i]) + Points.size(0) * nx) - 1];
                  }
                }

                y.set_size(1, predictedPoints.size(1));
                for (nx = 0; nx < c_loop_ub; nx++) {
                  y[nx] = predictedPoints[(validPredictions[b_i] +
                    predictedPoints.size(0) * nx) - 1];
                }

                bsxfun(c_this, y, diffs);
                dists.set_size(diffs.size(0));
                nx = diffs.size(0);
                for (loop_ub = 0; loop_ub < nx; loop_ub++) {
                  dists[loop_ub] = rt_hypotf_snf(diffs[loop_ub], diffs[loop_ub +
                    diffs.size(0)]);
                }

                loop_ub = indices.size(1) - 1;
                nx = 0;
                for (int c_i{0}; c_i <= loop_ub; c_i++) {
                  if (indices[c_i] > 0.0) {
                    nx++;
                  }
                }

                r1.set_size(1, nx);
                nx = 0;
                for (int c_i{0}; c_i <= loop_ub; c_i++) {
                  if (indices[c_i] > 0.0) {
                    r1[nx] = c_i + 1;
                    nx++;
                  }
                }

                y.set_size(1, r1.size(1));
                loop_ub = r1.size(1);
                for (nx = 0; nx < loop_ub; nx++) {
                  y[nx] = indices[r1[nx] - 1];
                }

                isMember(remIdx, y, distIdx);
                loop_ub = distIdx.size(1);
                for (int c_i{0}; c_i < loop_ub; c_i++) {
                  if (distIdx[c_i]) {
                    dists[c_i] = rtInfF;
                  }
                }

                ::coder::internal::minimum(dists, &ex, &nx);
                indices[validPredictions[b_i] - 1] = remIdx[nx - 1];
              }
            }
          }

          //
          // Arguments    : const double predictedPoint[2]
          //                double radius
          //                const double curve_data[]
          //                const int curve_size[2]
          //                const double coordsToUse[2]
          //                const ::coder::array<double, 2U> &removedIdx
          //                ::coder::array<double, 2U> &idx
          // Return Type  : void
          //
          void Checkerboard::findClosestOnCurve(const double predictedPoint[2],
            double radius, const double curve_data[], const int curve_size[2],
            const double coordsToUse[2], const ::coder::array<double, 2U>
            &removedIdx, ::coder::array<double, 2U> &idx) const
          {
            array<double, 2U> dataPts;
            array<double, 2U> firstCoord;
            array<double, 2U> remIdx;
            array<double, 2U> y;
            array<double, 1U> b_this;
            array<double, 1U> dist;
            array<float, 2U> b_a;
            array<float, 2U> currPt;
            array<float, 2U> diffs;
            array<float, 2U> queryPts;
            array<float, 1U> dists;
            array<float, 1U> varargin_1;
            array<int, 1U> ii;
            array<int, 1U> r;
            array<boolean_T, 1U> s;
            int acoef;
            int b_acoef;
            int b_k;
            int b_loop_ub;
            int c_k;
            int c_loop_ub;
            int i;
            int i1;
            int loop_ub;
            int ncols;
            int outsize_idx_1;
            acoef = Points.size(0);
            if (acoef < 1) {
              y.set_size(1, 0);
            } else {
              y.set_size(1, acoef);
              loop_ub = acoef - 1;
              if (static_cast<int>(acoef < 3200)) {
                for (i = 0; i <= loop_ub; i++) {
                  y[i] = static_cast<double>(i) + 1.0;
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (i = 0; i <= loop_ub; i++) {
                  y[i] = static_cast<double>(i) + 1.0;
                }
              }
            }

            b_acoef = BoardIdx.size(0) * BoardIdx.size(1);
            b_this = BoardIdx.reshape(b_acoef);
            do_vectors(y, b_this, remIdx, ii, &acoef);
            y.set_size(1, remIdx.size(1));
            loop_ub = remIdx.size(0) * remIdx.size(1) - 1;
            for (i1 = 0; i1 <= loop_ub; i1++) {
              y[i1] = remIdx[i1];
            }

            do_vectors(y, removedIdx, remIdx, ii, &acoef);
            diffs.set_size(remIdx.size(1), 2);
            if (remIdx.size(1) != 0) {
              b_acoef = (Points.size(1) != 1);
              acoef = (remIdx.size(1) != 1);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(b_k,i,b_loop_ub)

              for (int k = 0; k < 2; k++) {
                b_loop_ub = b_acoef * k;
                i = diffs.size(0) - 1;
                for (b_k = 0; b_k <= i; b_k++) {
                  diffs[b_k + diffs.size(0) * k] = Points[(static_cast<int>
                    (remIdx[acoef * b_k]) + Points.size(0) * b_loop_ub) - 1] -
                    static_cast<float>(predictedPoint[k]);
                }
              }
            }

            dists.set_size(diffs.size(0));
            b_acoef = diffs.size(0);
            if (static_cast<int>(diffs.size(0) < 3200)) {
              for (int k{0}; k < b_acoef; k++) {
                dists[k] = rt_hypotf_snf(diffs[k], diffs[k + diffs.size(0)]);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int k = 0; k < b_acoef; k++) {
                dists[k] = rt_hypotf_snf(diffs[k], diffs[k + diffs.size(0)]);
              }
            }

            s.set_size(dists.size(0));
            loop_ub = dists.size(0);
            if (static_cast<int>(dists.size(0) < 3200)) {
              for (i = 0; i < loop_ub; i++) {
                s[i] = (dists[i] < radius);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (i = 0; i < loop_ub; i++) {
                s[i] = (dists[i] < radius);
              }
            }

            acoef = 0;
            i1 = s.size(0);
            for (c_k = 0; c_k < i1; c_k++) {
              if (s[c_k]) {
                acoef++;
              }
            }

            if (acoef > 1) {
              double a;
              double a_tmp;
              int ntilerows;
              boolean_T exitg1;
              boolean_T p;
              a_tmp = predictedPoint[static_cast<int>(coordsToUse[0]) - 1];
              a = a_tmp - radius;
              a_tmp += radius;
              if (std::isnan(a) || std::isnan(a_tmp)) {
                firstCoord.set_size(1, 1);
                firstCoord[0] = rtNaN;
              } else if (a_tmp < a) {
                firstCoord.set_size(1, 0);
              } else if ((std::isinf(a) || std::isinf(a_tmp)) && (a == a_tmp)) {
                firstCoord.set_size(1, 1);
                firstCoord[0] = rtNaN;
              } else if (std::floor(a) == a) {
                loop_ub = static_cast<int>(a_tmp - a);
                firstCoord.set_size(1, loop_ub + 1);
                if (static_cast<int>(loop_ub + 1 < 3200)) {
                  for (i = 0; i <= loop_ub; i++) {
                    firstCoord[i] = a + static_cast<double>(i);
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i <= loop_ub; i++) {
                    firstCoord[i] = a + static_cast<double>(i);
                  }
                }
              } else {
                eml_float_colon(a, a_tmp, firstCoord);
              }

              p = true;
              c_k = 0;
              exitg1 = false;
              while ((!exitg1) && (c_k < 2)) {
                if (!(coordsToUse[c_k] == static_cast<double>(c_k) + 1.0)) {
                  p = false;
                  exitg1 = true;
                } else {
                  c_k++;
                }
              }

              if (p) {
                y.set_size(1, firstCoord.size(1));
                if (firstCoord.size(1) != 0) {
                  acoef = firstCoord.size(1);
                  y.set_size(1, firstCoord.size(1));
                  loop_ub = firstCoord.size(1);
                  if (static_cast<int>(firstCoord.size(1) < 3200)) {
                    for (i = 0; i < acoef; i++) {
                      y[i] = curve_data[0];
                    }
                  } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                    for (i = 0; i < loop_ub; i++) {
                      y[i] = curve_data[0];
                    }
                  }

                  i1 = curve_size[1];
                  for (c_k = 0; c_k <= i1 - 2; c_k++) {
                    if (firstCoord.size(1) == y.size(1)) {
                      a_tmp = curve_data[c_k + 1];
                      loop_ub = firstCoord.size(1) - 1;
                      y.set_size(1, firstCoord.size(1));
                      for (b_acoef = 0; b_acoef <= loop_ub; b_acoef++) {
                        y[b_acoef] = firstCoord[b_acoef] * y[b_acoef] + a_tmp;
                      }
                    } else {
                      binary_expand_op(y, firstCoord, curve_data, c_k);
                    }
                  }
                }

                dataPts.set_size(firstCoord.size(1), 2);
                loop_ub = firstCoord.size(1);
                if (static_cast<int>(firstCoord.size(1) < 3200)) {
                  for (i = 0; i < loop_ub; i++) {
                    dataPts[i] = firstCoord[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i < loop_ub; i++) {
                    dataPts[i] = firstCoord[i];
                  }
                }

                loop_ub = y.size(1);
                if (static_cast<int>(y.size(1) < 3200)) {
                  for (i = 0; i < loop_ub; i++) {
                    dataPts[i + dataPts.size(0)] = y[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i < loop_ub; i++) {
                    dataPts[i + dataPts.size(0)] = y[i];
                  }
                }
              } else {
                y.set_size(1, firstCoord.size(1));
                if (firstCoord.size(1) != 0) {
                  acoef = firstCoord.size(1);
                  y.set_size(1, firstCoord.size(1));
                  loop_ub = firstCoord.size(1);
                  if (static_cast<int>(firstCoord.size(1) < 3200)) {
                    for (i = 0; i < acoef; i++) {
                      y[i] = curve_data[0];
                    }
                  } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                    for (i = 0; i < loop_ub; i++) {
                      y[i] = curve_data[0];
                    }
                  }

                  i1 = curve_size[1];
                  for (c_k = 0; c_k <= i1 - 2; c_k++) {
                    if (firstCoord.size(1) == y.size(1)) {
                      a_tmp = curve_data[c_k + 1];
                      loop_ub = firstCoord.size(1) - 1;
                      y.set_size(1, firstCoord.size(1));
                      for (b_acoef = 0; b_acoef <= loop_ub; b_acoef++) {
                        y[b_acoef] = firstCoord[b_acoef] * y[b_acoef] + a_tmp;
                      }
                    } else {
                      binary_expand_op(y, firstCoord, curve_data, c_k);
                    }
                  }
                }

                dataPts.set_size(y.size(1), 2);
                loop_ub = y.size(1);
                if (static_cast<int>(y.size(1) < 3200)) {
                  for (i = 0; i < loop_ub; i++) {
                    dataPts[i] = y[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i < loop_ub; i++) {
                    dataPts[i] = y[i];
                  }
                }

                loop_ub = firstCoord.size(1);
                if (static_cast<int>(firstCoord.size(1) < 3200)) {
                  for (i = 0; i < loop_ub; i++) {
                    dataPts[i + dataPts.size(0)] = firstCoord[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i < loop_ub; i++) {
                    dataPts[i + dataPts.size(0)] = firstCoord[i];
                  }
                }
              }

              b_acoef = dists.size(0) - 1;
              acoef = 0;
              for (c_k = 0; c_k <= b_acoef; c_k++) {
                if (dists[c_k] < radius) {
                  acoef++;
                }
              }

              r.set_size(acoef);
              acoef = 0;
              for (c_k = 0; c_k <= b_acoef; c_k++) {
                if (dists[c_k] < radius) {
                  r[acoef] = c_k + 1;
                  acoef++;
                }
              }

              queryPts.set_size(r.size(0), Points.size(1));
              loop_ub = Points.size(1);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(b_k,b_loop_ub)

              for (i = 0; i < loop_ub; i++) {
                b_loop_ub = r.size(0);
                for (b_k = 0; b_k < b_loop_ub; b_k++) {
                  queryPts[b_k + queryPts.size(0) * i] = Points[(static_cast<int>
                    (remIdx[r[b_k] - 1]) + Points.size(0) * i) - 1];
                }
              }

              dist.set_size(r.size(0));
              i1 = r.size(0);
              if (r.size(0) - 1 >= 0) {
                c_loop_ub = Points.size(1);
                outsize_idx_1 = Points.size(1);
                ncols = Points.size(1);
              }

              for (c_k = 0; c_k < i1; c_k++) {
                b_a.set_size(1, Points.size(1));
                for (b_acoef = 0; b_acoef < c_loop_ub; b_acoef++) {
                  b_a[b_acoef] = queryPts[c_k + queryPts.size(0) * b_acoef];
                }

                currPt.set_size(dataPts.size(0), outsize_idx_1);
                ntilerows = dataPts.size(0);
                for (acoef = 0; acoef < ncols; acoef++) {
                  b_acoef = acoef * ntilerows;
                  for (loop_ub = 0; loop_ub < ntilerows; loop_ub++) {
                    currPt[b_acoef + loop_ub] = b_a[acoef];
                  }
                }

                if ((dataPts.size(0) == currPt.size(0)) && (currPt.size(1) == 2))
                {
                  diffs.set_size(dataPts.size(0), 2);
                  loop_ub = dataPts.size(0) * 2;
                  for (b_acoef = 0; b_acoef < loop_ub; b_acoef++) {
                    float b_varargin_1;
                    b_varargin_1 = static_cast<float>(dataPts[b_acoef]) -
                      currPt[b_acoef];
                    diffs[b_acoef] = b_varargin_1 * b_varargin_1;
                  }
                } else {
                  binary_expand_op(diffs, dataPts, currPt);
                }

                if (diffs.size(0) == 0) {
                  varargin_1.set_size(0);
                } else {
                  acoef = diffs.size(0);
                  varargin_1.set_size(diffs.size(0));
                  for (b_acoef = 0; b_acoef < acoef; b_acoef++) {
                    varargin_1[b_acoef] = diffs[b_acoef];
                  }

                  for (b_acoef = 0; b_acoef < acoef; b_acoef++) {
                    varargin_1[b_acoef] = varargin_1[b_acoef] + diffs[acoef +
                      b_acoef];
                  }
                }

                dist[c_k] = std::sqrt(::coder::internal::minimum(varargin_1));
              }

              b_acoef = dist.size(0);
              if (dist.size(0) <= 2) {
                if (dist.size(0) == 1) {
                  ntilerows = 1;
                } else if ((dist[0] > dist[dist.size(0) - 1]) || (std::isnan
                            (dist[0]) && (!std::isnan(dist[dist.size(0) - 1]))))
                {
                  ntilerows = dist.size(0);
                } else {
                  ntilerows = 1;
                }
              } else {
                if (!std::isnan(dist[0])) {
                  ntilerows = 1;
                } else {
                  ntilerows = 0;
                  c_k = 2;
                  exitg1 = false;
                  while ((!exitg1) && (c_k <= b_acoef)) {
                    if (!std::isnan(dist[c_k - 1])) {
                      ntilerows = c_k;
                      exitg1 = true;
                    } else {
                      c_k++;
                    }
                  }
                }

                if (ntilerows == 0) {
                  ntilerows = 1;
                } else {
                  a_tmp = dist[ntilerows - 1];
                  i1 = ntilerows + 1;
                  for (c_k = i1; c_k <= b_acoef; c_k++) {
                    a = dist[c_k - 1];
                    if (a_tmp > a) {
                      a_tmp = a;
                      ntilerows = c_k;
                    }
                  }
                }
              }

              s.set_size(dists.size(0));
              loop_ub = dists.size(0);
              if (static_cast<int>(dists.size(0) < 3200)) {
                for (i = 0; i < loop_ub; i++) {
                  s[i] = (dists[i] < radius);
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (i = 0; i < loop_ub; i++) {
                  s[i] = (dists[i] < radius);
                }
              }

              b_acoef = s.size(0);
              if (ntilerows <= b_acoef) {
                b_acoef = ntilerows;
              }

              ntilerows = 0;
              ii.set_size(b_acoef);
              acoef = 0;
              exitg1 = false;
              while ((!exitg1) && (acoef <= s.size(0) - 1)) {
                if (s[acoef]) {
                  ntilerows++;
                  ii[ntilerows - 1] = acoef + 1;
                  if (ntilerows >= b_acoef) {
                    exitg1 = true;
                  } else {
                    acoef++;
                  }
                } else {
                  acoef++;
                }
              }

              if (b_acoef == 1) {
                if (ntilerows == 0) {
                  ii.set_size(0);
                }
              } else {
                if (ntilerows < 1) {
                  ntilerows = 0;
                }

                ii.set_size(ntilerows);
              }

              dist.set_size(ii.size(0));
              loop_ub = ii.size(0);
              if (static_cast<int>(ii.size(0) < 3200)) {
                for (i = 0; i < loop_ub; i++) {
                  dist[i] = ii[i];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (i = 0; i < loop_ub; i++) {
                  dist[i] = ii[i];
                }
              }

              idx.set_size(1, 1);
              idx[0] = remIdx[static_cast<int>(dist[dist.size(0) - 1]) - 1];
            } else {
              s.set_size(dists.size(0));
              loop_ub = dists.size(0);
              if (static_cast<int>(dists.size(0) < 3200)) {
                for (i = 0; i < loop_ub; i++) {
                  s[i] = (dists[i] < radius);
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (i = 0; i < loop_ub; i++) {
                  s[i] = (dists[i] < radius);
                }
              }

              acoef = 0;
              i1 = s.size(0);
              for (c_k = 0; c_k < i1; c_k++) {
                if (s[c_k]) {
                  acoef++;
                }
              }

              if (acoef == 1) {
                b_acoef = dists.size(0) - 1;
                acoef = 0;
                for (c_k = 0; c_k <= b_acoef; c_k++) {
                  if (dists[c_k] < radius) {
                    acoef++;
                  }
                }

                r.set_size(acoef);
                acoef = 0;
                for (c_k = 0; c_k <= b_acoef; c_k++) {
                  if (dists[c_k] < radius) {
                    r[acoef] = c_k + 1;
                    acoef++;
                  }
                }

                idx.set_size(1, r.size(0));
                loop_ub = r.size(0);
                if (static_cast<int>(r.size(0) < 3200)) {
                  for (i = 0; i < loop_ub; i++) {
                    idx[i] = remIdx[r[i] - 1];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i < loop_ub; i++) {
                    idx[i] = remIdx[r[i] - 1];
                  }
                }
              } else {
                idx.set_size(1, 0);
              }
            }
          }

          //
          // Arguments    : const ::coder::array<double, 2U> &idx
          //                double coordsToUse[2]
          // Return Type  : void
          //
          void Checkerboard::findIndependentVar(const ::coder::array<double, 2U>
            &idx, double coordsToUse[2]) const
          {
            array<double, 2U> b_x;
            array<double, 2U> r4;
            array<double, 2U> x;
            array<int, 2U> r2;
            array<int, 2U> r3;
            array<int, 2U> r5;
            array<int, 2U> r6;
            array<boolean_T, 2U> r;
            array<boolean_T, 2U> r1;
            int idx_tmp;
            int loop_ub;
            idx_tmp = static_cast<int>(idx[0]);
            r.set_size(1, BoardIdx.size(1));
            loop_ub = BoardIdx.size(1);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                r[i] = (BoardIdx[(idx_tmp + BoardIdx.size(0) * i) - 1] > 0.0);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                r[i] = (BoardIdx[(idx_tmp + BoardIdx.size(0) * i) - 1] > 0.0);
              }
            }

            idx_tmp = static_cast<int>(idx[1]);
            r1.set_size(1, BoardIdx.size(1));
            loop_ub = BoardIdx.size(1);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                r1[i] = (BoardIdx[(idx_tmp + BoardIdx.size(0) * i) - 1] > 0.0);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                r1[i] = (BoardIdx[(idx_tmp + BoardIdx.size(0) * i) - 1] > 0.0);
              }
            }

            loop_ub = r.size(1) - 1;
            idx_tmp = 0;
            for (int b_i{0}; b_i <= loop_ub; b_i++) {
              if (r[b_i] && r1[b_i]) {
                idx_tmp++;
              }
            }

            r2.set_size(1, idx_tmp);
            idx_tmp = 0;
            for (int b_i{0}; b_i <= loop_ub; b_i++) {
              if (r[b_i] && r1[b_i]) {
                r2[idx_tmp] = b_i + 1;
                idx_tmp++;
              }
            }

            idx_tmp = static_cast<int>(idx[1]);
            x.set_size(1, r2.size(1));
            loop_ub = r2.size(1);
            if (static_cast<int>(r2.size(1) < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                x[i] = BoardCoords[(idx_tmp + BoardCoords.size(0) * (r2[i] - 1))
                  - 1];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                x[i] = BoardCoords[(idx_tmp + BoardCoords.size(0) * (r2[i] - 1))
                  - 1];
              }
            }

            loop_ub = r.size(1) - 1;
            idx_tmp = 0;
            for (int b_i{0}; b_i <= loop_ub; b_i++) {
              if (r[b_i] && r1[b_i]) {
                idx_tmp++;
              }
            }

            r3.set_size(1, idx_tmp);
            idx_tmp = 0;
            for (int b_i{0}; b_i <= loop_ub; b_i++) {
              if (r[b_i] && r1[b_i]) {
                r3[idx_tmp] = b_i + 1;
                idx_tmp++;
              }
            }

            idx_tmp = static_cast<int>(idx[0]);
            r4.set_size(1, r3.size(1));
            loop_ub = r3.size(1);
            if (static_cast<int>(r3.size(1) < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                r4[i] = BoardCoords[(idx_tmp + BoardCoords.size(0) * (r3[i] - 1))
                  - 1];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                r4[i] = BoardCoords[(idx_tmp + BoardCoords.size(0) * (r3[i] - 1))
                  - 1];
              }
            }

            x.set_size(1, x.size(1));
            idx_tmp = x.size(1) - 1;
            loop_ub = x.size(1) - 1;
            if (static_cast<int>(x.size(1) < 3200)) {
              for (int i{0}; i <= idx_tmp; i++) {
                x[i] = x[i] - r4[i];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i <= loop_ub; i++) {
                x[i] = x[i] - r4[i];
              }
            }

            loop_ub = r.size(1) - 1;
            idx_tmp = 0;
            for (int b_i{0}; b_i <= loop_ub; b_i++) {
              if (r[b_i] && r1[b_i]) {
                idx_tmp++;
              }
            }

            r5.set_size(1, idx_tmp);
            idx_tmp = 0;
            for (int b_i{0}; b_i <= loop_ub; b_i++) {
              if (r[b_i] && r1[b_i]) {
                r5[idx_tmp] = b_i + 1;
                idx_tmp++;
              }
            }

            idx_tmp = static_cast<int>(idx[1]);
            b_x.set_size(1, r5.size(1));
            loop_ub = r5.size(1);
            if (static_cast<int>(r5.size(1) < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                b_x[i] = BoardCoords[((idx_tmp + BoardCoords.size(0) * (r5[i] -
                  1)) + BoardCoords.size(0) * BoardCoords.size(1)) - 1];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                b_x[i] = BoardCoords[((idx_tmp + BoardCoords.size(0) * (r5[i] -
                  1)) + BoardCoords.size(0) * BoardCoords.size(1)) - 1];
              }
            }

            loop_ub = r.size(1) - 1;
            idx_tmp = 0;
            for (int b_i{0}; b_i <= loop_ub; b_i++) {
              if (r[b_i] && r1[b_i]) {
                idx_tmp++;
              }
            }

            r6.set_size(1, idx_tmp);
            idx_tmp = 0;
            for (int b_i{0}; b_i <= loop_ub; b_i++) {
              if (r[b_i] && r1[b_i]) {
                r6[idx_tmp] = b_i + 1;
                idx_tmp++;
              }
            }

            idx_tmp = static_cast<int>(idx[0]);
            r4.set_size(1, r6.size(1));
            loop_ub = r6.size(1);
            if (static_cast<int>(r6.size(1) < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                r4[i] = BoardCoords[((idx_tmp + BoardCoords.size(0) * (r6[i] - 1))
                                     + BoardCoords.size(0) * BoardCoords.size(1))
                  - 1];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                r4[i] = BoardCoords[((idx_tmp + BoardCoords.size(0) * (r6[i] - 1))
                                     + BoardCoords.size(0) * BoardCoords.size(1))
                  - 1];
              }
            }

            b_x.set_size(1, b_x.size(1));
            idx_tmp = b_x.size(1) - 1;
            loop_ub = b_x.size(1) - 1;
            if (static_cast<int>(b_x.size(1) < 3200)) {
              for (int i{0}; i <= idx_tmp; i++) {
                b_x[i] = b_x[i] - r4[i];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i <= loop_ub; i++) {
                b_x[i] = b_x[i] - r4[i];
              }
            }

            if (std::abs(combineVectorElements(x) / static_cast<double>(x.size(1)))
                > std::abs(combineVectorElements(b_x) / static_cast<double>
                           (b_x.size(1)))) {
              coordsToUse[0] = 1.0;
              coordsToUse[1] = 2.0;
            } else {
              coordsToUse[0] = 2.0;
              coordsToUse[1] = 1.0;
            }
          }

          //
          // Arguments    : double coordsToUse[2]
          // Return Type  : void
          //
          void Checkerboard::findIndependentVar(double coordsToUse[2]) const
          {
            array<double, 2U> b_x;
            array<double, 2U> r4;
            array<double, 2U> x;
            array<int, 2U> r2;
            array<int, 2U> r3;
            array<int, 2U> r5;
            array<int, 2U> r6;
            array<boolean_T, 2U> r;
            array<boolean_T, 2U> r1;
            int end;
            int loop_ub;
            r.set_size(1, BoardIdx.size(1));
            loop_ub = BoardIdx.size(1);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                r[i] = (BoardIdx[BoardIdx.size(0) * i] > 0.0);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                r[i] = (BoardIdx[BoardIdx.size(0) * i] > 0.0);
              }
            }

            r1.set_size(1, BoardIdx.size(1));
            loop_ub = BoardIdx.size(1);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                r1[i] = (BoardIdx[BoardIdx.size(0) * i + 1] > 0.0);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                r1[i] = (BoardIdx[BoardIdx.size(0) * i + 1] > 0.0);
              }
            }

            end = r.size(1) - 1;
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (r[b_i] && r1[b_i]) {
                loop_ub++;
              }
            }

            r2.set_size(1, loop_ub);
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (r[b_i] && r1[b_i]) {
                r2[loop_ub] = b_i + 1;
                loop_ub++;
              }
            }

            x.set_size(1, r2.size(1));
            loop_ub = r2.size(1);
            if (static_cast<int>(r2.size(1) < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                x[i] = BoardCoords[BoardCoords.size(0) * (r2[i] - 1) + 1];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                x[i] = BoardCoords[BoardCoords.size(0) * (r2[i] - 1) + 1];
              }
            }

            end = r.size(1) - 1;
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (r[b_i] && r1[b_i]) {
                loop_ub++;
              }
            }

            r3.set_size(1, loop_ub);
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (r[b_i] && r1[b_i]) {
                r3[loop_ub] = b_i + 1;
                loop_ub++;
              }
            }

            r4.set_size(1, r3.size(1));
            loop_ub = r3.size(1);
            if (static_cast<int>(r3.size(1) < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                r4[i] = BoardCoords[BoardCoords.size(0) * (r3[i] - 1)];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                r4[i] = BoardCoords[BoardCoords.size(0) * (r3[i] - 1)];
              }
            }

            x.set_size(1, x.size(1));
            end = x.size(1) - 1;
            loop_ub = x.size(1) - 1;
            if (static_cast<int>(x.size(1) < 3200)) {
              for (int i{0}; i <= end; i++) {
                x[i] = x[i] - r4[i];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i <= loop_ub; i++) {
                x[i] = x[i] - r4[i];
              }
            }

            end = r.size(1) - 1;
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (r[b_i] && r1[b_i]) {
                loop_ub++;
              }
            }

            r5.set_size(1, loop_ub);
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (r[b_i] && r1[b_i]) {
                r5[loop_ub] = b_i + 1;
                loop_ub++;
              }
            }

            b_x.set_size(1, r5.size(1));
            loop_ub = r5.size(1);
            if (static_cast<int>(r5.size(1) < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                b_x[i] = BoardCoords[(BoardCoords.size(0) * (r5[i] - 1) +
                                      BoardCoords.size(0) * BoardCoords.size(1))
                  + 1];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                b_x[i] = BoardCoords[(BoardCoords.size(0) * (r5[i] - 1) +
                                      BoardCoords.size(0) * BoardCoords.size(1))
                  + 1];
              }
            }

            end = r.size(1) - 1;
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (r[b_i] && r1[b_i]) {
                loop_ub++;
              }
            }

            r6.set_size(1, loop_ub);
            loop_ub = 0;
            for (int b_i{0}; b_i <= end; b_i++) {
              if (r[b_i] && r1[b_i]) {
                r6[loop_ub] = b_i + 1;
                loop_ub++;
              }
            }

            r4.set_size(1, r6.size(1));
            loop_ub = r6.size(1);
            if (static_cast<int>(r6.size(1) < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                r4[i] = BoardCoords[BoardCoords.size(0) * (r6[i] - 1) +
                  BoardCoords.size(0) * BoardCoords.size(1)];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                r4[i] = BoardCoords[BoardCoords.size(0) * (r6[i] - 1) +
                  BoardCoords.size(0) * BoardCoords.size(1)];
              }
            }

            b_x.set_size(1, b_x.size(1));
            end = b_x.size(1) - 1;
            loop_ub = b_x.size(1) - 1;
            if (static_cast<int>(b_x.size(1) < 3200)) {
              for (int i{0}; i <= end; i++) {
                b_x[i] = b_x[i] - r4[i];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i <= loop_ub; i++) {
                b_x[i] = b_x[i] - r4[i];
              }
            }

            if (std::abs(combineVectorElements(x) / static_cast<double>(x.size(1)))
                > std::abs(combineVectorElements(b_x) / static_cast<double>
                           (b_x.size(1)))) {
              coordsToUse[0] = 1.0;
              coordsToUse[1] = 2.0;
            } else {
              coordsToUse[0] = 2.0;
              coordsToUse[1] = 1.0;
            }
          }

          //
          // Arguments    : const ::coder::array<float, 2U> &pointVectors
          //                const ::coder::array<float, 1U> &euclideanDists
          //                const ::coder::array<float, 2U> &v
          // Return Type  : double
          //
          double Checkerboard::findNeighbor(const ::coder::array<float, 2U>
            &pointVectors, const ::coder::array<float, 1U> &euclideanDists,
            const ::coder::array<float, 2U> &v) const
          {
            array<float, 1U> angleCosines;
            array<float, 1U> dists;
            array<int, 1U> r1;
            array<int, 1U> r2;
            array<boolean_T, 2U> r;
            double neighborIdx;
            float b;
            int aoffset;
            int inner;
            int mc;
            mc = pointVectors.size(0) - 1;
            inner = pointVectors.size(1);
            angleCosines.set_size(pointVectors.size(0));
            if (static_cast<int>(pointVectors.size(0) < 3200)) {
              for (int i{0}; i <= mc; i++) {
                angleCosines[i] = 0.0F;
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i <= mc; i++) {
                angleCosines[i] = 0.0F;
              }
            }

            for (int k{0}; k < inner; k++) {
              aoffset = k * pointVectors.size(0);
              for (int b_i{0}; b_i <= mc; b_i++) {
                angleCosines[b_i] = angleCosines[b_i] + pointVectors[aoffset +
                  b_i] * v[k];
              }
            }

            b = rt_hypotf_snf(v[0], v[1]);
            if (angleCosines.size(0) == euclideanDists.size(0)) {
              aoffset = angleCosines.size(0);
              if (static_cast<int>(angleCosines.size(0) < 3200)) {
                for (int i{0}; i < aoffset; i++) {
                  angleCosines[i] = angleCosines[i] / (euclideanDists[i] * b);
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < aoffset; i++) {
                  angleCosines[i] = angleCosines[i] / (euclideanDists[i] * b);
                }
              }
            } else {
              binary_expand_op(angleCosines, euclideanDists, b);
            }

            if (euclideanDists.size(0) == 1) {
              inner = angleCosines.size(0);
            } else {
              inner = euclideanDists.size(0);
            }

            if ((euclideanDists.size(0) == angleCosines.size(0)) &&
                (euclideanDists.size(0) == inner)) {
              dists.set_size(euclideanDists.size(0));
              aoffset = euclideanDists.size(0);
              if (static_cast<int>(euclideanDists.size(0) < 3200)) {
                for (int i{0}; i < aoffset; i++) {
                  dists[i] = euclideanDists[i] + 1.5F * euclideanDists[i] *
                    (1.0F - angleCosines[i]);
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < aoffset; i++) {
                  dists[i] = euclideanDists[i] + 1.5F * euclideanDists[i] *
                    (1.0F - angleCosines[i]);
                }
              }
            } else {
              binary_expand_op(dists, euclideanDists, angleCosines);
            }

            r.set_size(BoardIdx.size(0), BoardIdx.size(1));
            aoffset = BoardIdx.size(0) * BoardIdx.size(1);
            if (static_cast<int>(aoffset < 3200)) {
              for (int i{0}; i < aoffset; i++) {
                r[i] = (BoardIdx[i] > 0.0);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < aoffset; i++) {
                r[i] = (BoardIdx[i] > 0.0);
              }
            }

            mc = r.size(0) * r.size(1) - 1;
            aoffset = 0;
            for (int b_i{0}; b_i <= mc; b_i++) {
              if (r[b_i]) {
                aoffset++;
              }
            }

            r1.set_size(aoffset);
            aoffset = 0;
            for (int b_i{0}; b_i <= mc; b_i++) {
              if (r[b_i]) {
                r1[aoffset] = b_i + 1;
                aoffset++;
              }
            }

            r2.set_size(r1.size(0));
            aoffset = r1.size(0);
            if (static_cast<int>(r1.size(0) < 3200)) {
              for (int i{0}; i < aoffset; i++) {
                r2[i] = static_cast<int>(BoardIdx[r1[i] - 1]);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < aoffset; i++) {
                r2[i] = static_cast<int>(BoardIdx[r1[i] - 1]);
              }
            }

            aoffset = r2.size(0);
            for (mc = 0; mc < aoffset; mc++) {
              dists[r2[mc] - 1] = rtInfF;
            }

            mc = angleCosines.size(0);
            if (static_cast<int>(angleCosines.size(0) < 3200)) {
              for (int i{0}; i < mc; i++) {
                if (angleCosines[i] < 0.0F) {
                  dists[i] = rtInfF;
                }
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < mc; i++) {
                if (angleCosines[i] < 0.0F) {
                  dists[i] = rtInfF;
                }
              }
            }

            ::coder::internal::minimum(dists, &b, &aoffset);
            neighborIdx = aoffset;
            if (std::isinf(b)) {
              neighborIdx = -1.0;
            }

            return neighborIdx;
          }

          //
          // Arguments    : const ::coder::array<float, 2U> &pointVectors
          //                const ::coder::array<float, 1U> &euclideanDists
          //                const float v[2]
          // Return Type  : double
          //
          double Checkerboard::findNeighbor(const ::coder::array<float, 2U>
            &pointVectors, const ::coder::array<float, 1U> &euclideanDists,
            const float v[2]) const
          {
            array<float, 1U> angleCosines;
            array<float, 1U> dists;
            array<int, 1U> r1;
            array<int, 1U> r2;
            array<boolean_T, 2U> r;
            double neighborIdx;
            float b;
            int aoffset;
            int inner;
            int mc;
            mc = pointVectors.size(0) - 1;
            inner = pointVectors.size(1);
            angleCosines.set_size(pointVectors.size(0));
            if (static_cast<int>(pointVectors.size(0) < 3200)) {
              for (int i{0}; i <= mc; i++) {
                angleCosines[i] = 0.0F;
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i <= mc; i++) {
                angleCosines[i] = 0.0F;
              }
            }

            for (int k{0}; k < inner; k++) {
              aoffset = k * pointVectors.size(0);
              for (int b_i{0}; b_i <= mc; b_i++) {
                angleCosines[b_i] = angleCosines[b_i] + pointVectors[aoffset +
                  b_i] * v[k];
              }
            }

            b = rt_hypotf_snf(v[0], v[1]);
            if (angleCosines.size(0) == euclideanDists.size(0)) {
              aoffset = angleCosines.size(0);
              if (static_cast<int>(angleCosines.size(0) < 3200)) {
                for (int i{0}; i < aoffset; i++) {
                  angleCosines[i] = angleCosines[i] / (euclideanDists[i] * b);
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < aoffset; i++) {
                  angleCosines[i] = angleCosines[i] / (euclideanDists[i] * b);
                }
              }
            } else {
              binary_expand_op(angleCosines, euclideanDists, b);
            }

            if (euclideanDists.size(0) == 1) {
              inner = angleCosines.size(0);
            } else {
              inner = euclideanDists.size(0);
            }

            if ((euclideanDists.size(0) == angleCosines.size(0)) &&
                (euclideanDists.size(0) == inner)) {
              dists.set_size(euclideanDists.size(0));
              aoffset = euclideanDists.size(0);
              if (static_cast<int>(euclideanDists.size(0) < 3200)) {
                for (int i{0}; i < aoffset; i++) {
                  dists[i] = euclideanDists[i] + 1.5F * euclideanDists[i] *
                    (1.0F - angleCosines[i]);
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (int i = 0; i < aoffset; i++) {
                  dists[i] = euclideanDists[i] + 1.5F * euclideanDists[i] *
                    (1.0F - angleCosines[i]);
                }
              }
            } else {
              binary_expand_op(dists, euclideanDists, angleCosines);
            }

            r.set_size(BoardIdx.size(0), BoardIdx.size(1));
            aoffset = BoardIdx.size(0) * BoardIdx.size(1);
            if (static_cast<int>(aoffset < 3200)) {
              for (int i{0}; i < aoffset; i++) {
                r[i] = (BoardIdx[i] > 0.0);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < aoffset; i++) {
                r[i] = (BoardIdx[i] > 0.0);
              }
            }

            mc = r.size(0) * r.size(1) - 1;
            aoffset = 0;
            for (int b_i{0}; b_i <= mc; b_i++) {
              if (r[b_i]) {
                aoffset++;
              }
            }

            r1.set_size(aoffset);
            aoffset = 0;
            for (int b_i{0}; b_i <= mc; b_i++) {
              if (r[b_i]) {
                r1[aoffset] = b_i + 1;
                aoffset++;
              }
            }

            r2.set_size(r1.size(0));
            aoffset = r1.size(0);
            if (static_cast<int>(r1.size(0) < 3200)) {
              for (int i{0}; i < aoffset; i++) {
                r2[i] = static_cast<int>(BoardIdx[r1[i] - 1]);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < aoffset; i++) {
                r2[i] = static_cast<int>(BoardIdx[r1[i] - 1]);
              }
            }

            aoffset = r2.size(0);
            for (mc = 0; mc < aoffset; mc++) {
              dists[r2[mc] - 1] = rtInfF;
            }

            mc = angleCosines.size(0);
            if (static_cast<int>(angleCosines.size(0) < 3200)) {
              for (int i{0}; i < mc; i++) {
                if (angleCosines[i] < 0.0F) {
                  dists[i] = rtInfF;
                }
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < mc; i++) {
                if (angleCosines[i] < 0.0F) {
                  dists[i] = rtInfF;
                }
              }
            }

            ::coder::internal::minimum(dists, &b, &aoffset);
            neighborIdx = aoffset;
            if (std::isinf(b)) {
              neighborIdx = -1.0;
            }

            return neighborIdx;
          }

          //
          // Arguments    : const ::coder::array<double, 2U> &idx
          //                const ::coder::array<double, 1U> &validIdx
          //                double currIdx
          //                const double coordsToUse[2]
          //                double *coordDist
          //                double *moveMultiplier
          //                double *firstValidIdx
          // Return Type  : void
          //
          void Checkerboard::findSearchParams(const ::coder::array<double, 2U>
            &idx, const ::coder::array<double, 1U> &validIdx, double currIdx,
            const double coordsToUse[2], double *coordDist, double
            *moveMultiplier, double *firstValidIdx) const
          {
            if (idx[0] == 1.0) {
              *moveMultiplier = validIdx[0];
              *firstValidIdx = validIdx[0];
              *coordDist = (BoardCoords[((static_cast<int>(validIdx[0]) +
                BoardCoords.size(0) * (static_cast<int>(currIdx) - 1)) +
                BoardCoords.size(0) * BoardCoords.size(1) * (static_cast<int>
                (coordsToUse[0]) - 1)) - 1] - BoardCoords[((static_cast<int>
                (validIdx[1]) + BoardCoords.size(0) * (static_cast<int>(currIdx)
                - 1)) + BoardCoords.size(0) * BoardCoords.size(1) * (
                static_cast<int>(coordsToUse[0]) - 1)) - 1]) / (validIdx[1] -
                validIdx[0]);
            } else {
              *moveMultiplier = (static_cast<double>(BoardCoords.size(0)) -
                                 validIdx[validIdx.size(0) - 1]) + 1.0;
              *firstValidIdx = validIdx[validIdx.size(0) - 1];
              *coordDist = (BoardCoords[((static_cast<int>
                (validIdx[validIdx.size(0) - 1]) + BoardCoords.size(0) * (
                static_cast<int>(currIdx) - 1)) + BoardCoords.size(0) *
                BoardCoords.size(1) * (static_cast<int>(coordsToUse[0]) - 1)) -
                            1] - BoardCoords[((static_cast<int>
                (validIdx[validIdx.size(0) - 2]) + BoardCoords.size(0) * (
                static_cast<int>(currIdx) - 1)) + BoardCoords.size(0) *
                BoardCoords.size(1) * (static_cast<int>(coordsToUse[0]) - 1)) -
                            1]) / (validIdx[validIdx.size(0) - 1] -
                                   validIdx[validIdx.size(0) - 2]);
            }
          }

          //
          // Arguments    : const ::coder::array<double, 2U> &idx
          //                const ::coder::array<double, 2U> &validIdx
          //                double currIdx
          //                const double coordsToUse[2]
          //                double *coordDist
          //                double *moveMultiplier
          //                double *firstValidIdx
          // Return Type  : void
          //
          void Checkerboard::findSearchParams(const ::coder::array<double, 2U>
            &idx, const ::coder::array<double, 2U> &validIdx, double currIdx,
            const double coordsToUse[2], double *coordDist, double
            *moveMultiplier, double *firstValidIdx) const
          {
            if (idx[0] == 1.0) {
              *moveMultiplier = validIdx[0];
              *firstValidIdx = validIdx[0];
              *coordDist = (BoardCoords[((static_cast<int>(currIdx) +
                BoardCoords.size(0) * (static_cast<int>(validIdx[0]) - 1)) +
                BoardCoords.size(0) * BoardCoords.size(1) * (static_cast<int>
                (coordsToUse[0]) - 1)) - 1] - BoardCoords[((static_cast<int>
                (currIdx) + BoardCoords.size(0) * (static_cast<int>(validIdx[1])
                - 1)) + BoardCoords.size(0) * BoardCoords.size(1) * (
                static_cast<int>(coordsToUse[0]) - 1)) - 1]) / (validIdx[1] -
                validIdx[0]);
            } else {
              *moveMultiplier = (static_cast<double>(BoardCoords.size(1)) -
                                 validIdx[validIdx.size(1) - 1]) + 1.0;
              *firstValidIdx = validIdx[validIdx.size(1) - 1];
              *coordDist = (BoardCoords[((static_cast<int>(currIdx) +
                BoardCoords.size(0) * (static_cast<int>(validIdx[validIdx.size(1)
                - 1]) - 1)) + BoardCoords.size(0) * BoardCoords.size(1) * (
                static_cast<int>(coordsToUse[0]) - 1)) - 1] - BoardCoords[((
                static_cast<int>(currIdx) + BoardCoords.size(0) * (static_cast<
                int>(validIdx[validIdx.size(1) - 2]) - 1)) + BoardCoords.size(0)
                * BoardCoords.size(1) * (static_cast<int>(coordsToUse[0]) - 1))
                            - 1]) / (validIdx[validIdx.size(1) - 1] -
                validIdx[validIdx.size(1) - 2]);
            }
          }

          //
          // Arguments    : const ::coder::array<double, 2U> &idx
          //                ::coder::array<double, 2U> &newIndices
          // Return Type  : void
          //
          void Checkerboard::fitPolynomialIndices(const ::coder::array<double,
            2U> &idx, ::coder::array<double, 2U> &newIndices) const
          {
            array<double, 2U> b_this;
            array<double, 2U> removedIdx;
            array<double, 2U> validIdx;
            array<int, 2U> r;
            array<int, 1U> currCurve_tmp;
            double currCurve_data[5];
            double coordsToUse[2];
            double coordDist;
            double moveDistMultiplier;
            double refCoordValue;
            int i1;
            int loop_ub;
            b_findIndependentVar(idx, coordsToUse);
            newIndices.set_size(1, BoardCoords.size(0));
            loop_ub = BoardCoords.size(0);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                newIndices[i] = 0.0;
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                newIndices[i] = 0.0;
              }
            }

            removedIdx.set_size(1, 0);
            i1 = BoardCoords.size(0);
            for (int j{0}; j < i1; j++) {
              int i2;
              validIdx.set_size(1, BoardCoords.size(1));
              loop_ub = BoardCoords.size(1);
              for (i2 = 0; i2 < loop_ub; i2++) {
                validIdx[i2] = BoardCoords[(j + BoardCoords.size(0) * i2) +
                  BoardCoords.size(0) * BoardCoords.size(1) * (static_cast<int>
                  (coordsToUse[0]) - 1)];
              }

              eml_find(validIdx, r);
              validIdx.set_size(1, r.size(1));
              loop_ub = r.size(1);
              for (i2 = 0; i2 < loop_ub; i2++) {
                validIdx[i2] = r[i2];
              }

              if (validIdx.size(1) >= 2) {
                double currCoord;
                double currRad;
                int currCurve_size[2];
                int n;
                boolean_T exitg1;
                findSearchParams(idx, validIdx, static_cast<double>(j) + 1.0,
                                 coordsToUse, &coordDist, &moveDistMultiplier,
                                 &refCoordValue);
                n = 0;
                i2 = validIdx.size(1);
                for (loop_ub = 0; loop_ub < i2; loop_ub++) {
                  if (static_cast<int>(validIdx[loop_ub]) != 0) {
                    n++;
                  }
                }

                currCurve_tmp.set_size(validIdx.size(1));
                loop_ub = validIdx.size(1);
                for (i2 = 0; i2 < loop_ub; i2++) {
                  currCurve_tmp[i2] = static_cast<int>(validIdx[i2]);
                }

                validIdx.set_size(1, currCurve_tmp.size(0));
                loop_ub = currCurve_tmp.size(0);
                for (i2 = 0; i2 < loop_ub; i2++) {
                  validIdx[i2] = BoardCoords[(j + BoardCoords.size(0) *
                    (currCurve_tmp[i2] - 1)) + BoardCoords.size(0) *
                    BoardCoords.size(1) * (static_cast<int>(coordsToUse[0]) - 1)];
                }

                b_this.set_size(1, currCurve_tmp.size(0));
                loop_ub = currCurve_tmp.size(0);
                for (i2 = 0; i2 < loop_ub; i2++) {
                  b_this[i2] = BoardCoords[(j + BoardCoords.size(0) *
                    (currCurve_tmp[i2] - 1)) + BoardCoords.size(0) *
                    BoardCoords.size(1) * (static_cast<int>(coordsToUse[1]) - 1)];
                }

                if (n > 5) {
                  i2 = 4;
                } else {
                  i2 = 2;
                }

                polyfit(validIdx, b_this, static_cast<double>(i2),
                        currCurve_data, currCurve_size);
                currRad = coordDist / 4.0;
                refCoordValue = BoardCoords[(j + BoardCoords.size(0) * (
                  static_cast<int>(refCoordValue) - 1)) + BoardCoords.size(0) *
                  BoardCoords.size(1) * (static_cast<int>(coordsToUse[0]) - 1)];
                currCoord = currRad + refCoordValue;
                exitg1 = false;
                while ((!exitg1) && (std::abs(currCoord - refCoordValue) <
                                     moveDistMultiplier * 1.5 * std::abs
                                     (coordDist))) {
                  double currPt[2];
                  boolean_T exitg2;
                  boolean_T p;
                  p = true;
                  loop_ub = 0;
                  exitg2 = false;
                  while ((!exitg2) && (loop_ub < 2)) {
                    if (!(coordsToUse[loop_ub] == static_cast<double>(loop_ub) +
                          1.0)) {
                      p = false;
                      exitg2 = true;
                    } else {
                      loop_ub++;
                    }
                  }

                  if (p) {
                    double y;
                    y = currCurve_data[0];
                    i2 = currCurve_size[1];
                    for (loop_ub = 0; loop_ub <= i2 - 2; loop_ub++) {
                      y = currCoord * y + currCurve_data[loop_ub + 1];
                    }

                    currPt[0] = currCoord;
                    currPt[1] = y;
                  } else {
                    double y;
                    y = currCurve_data[0];
                    i2 = currCurve_size[1];
                    for (loop_ub = 0; loop_ub <= i2 - 2; loop_ub++) {
                      y = currCoord * y + currCurve_data[loop_ub + 1];
                    }

                    currPt[0] = y;
                    currPt[1] = currCoord;
                  }

                  findClosestOnCurve(currPt, std::abs(currRad), currCurve_data,
                                     currCurve_size, coordsToUse, removedIdx,
                                     validIdx);
                  if (validIdx.size(1) != 0) {
                    newIndices[j] = validIdx[0];
                    i2 = removedIdx.size(1);
                    loop_ub = validIdx.size(1);
                    removedIdx.set_size(removedIdx.size(0), removedIdx.size(1) +
                                        validIdx.size(1));
                    for (n = 0; n < loop_ub; n++) {
                      removedIdx[i2 + n] = validIdx[n];
                    }

                    exitg1 = true;
                  } else {
                    currCoord += currRad;
                  }
                }
              }
            }
          }

          //
          // Arguments    : ::coder::array<double, 2U> &newIndices
          // Return Type  : void
          //
          void Checkerboard::fitPolynomialIndices(::coder::array<double, 2U>
            &newIndices) const
          {
            array<double, 2U> b_index;
            array<double, 2U> removedIdx;
            array<double, 1U> b_this;
            array<double, 1U> c_this;
            array<int, 1U> validIdx;
            double currCurve_data[5];
            double coordsToUse[2];
            int i1;
            int loop_ub;
            findIndependentVar(coordsToUse);
            newIndices.set_size(1, BoardCoords.size(1));
            loop_ub = BoardCoords.size(1);
            if (static_cast<int>(loop_ub < 3200)) {
              for (int i{0}; i < loop_ub; i++) {
                newIndices[i] = 0.0;
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int i = 0; i < loop_ub; i++) {
                newIndices[i] = 0.0;
              }
            }

            removedIdx.set_size(1, 0);
            i1 = BoardCoords.size(1);
            for (int j{0}; j < i1; j++) {
              int i2;
              b_this.set_size(BoardCoords.size(0));
              loop_ub = BoardCoords.size(0);
              for (i2 = 0; i2 < loop_ub; i2++) {
                b_this[i2] = BoardCoords[(i2 + BoardCoords.size(0) * j) +
                  BoardCoords.size(0) * BoardCoords.size(1) * (static_cast<int>
                  (coordsToUse[0]) - 1)];
              }

              eml_find(b_this, validIdx);
              if (validIdx.size(0) >= 2) {
                double coordDist;
                double currCoord;
                double currRad;
                double refCoordValue;
                int currCurve_size[2];
                int k;
                boolean_T exitg1;
                coordDist = (BoardCoords[((validIdx[0] + BoardCoords.size(0) * j)
                  + BoardCoords.size(0) * BoardCoords.size(1) * (static_cast<int>
                  (coordsToUse[0]) - 1)) - 1] - BoardCoords[((validIdx[1] +
                  BoardCoords.size(0) * j) + BoardCoords.size(0) *
                  BoardCoords.size(1) * (static_cast<int>(coordsToUse[0]) - 1))
                             - 1]) / (static_cast<double>(validIdx[1]) -
                                      static_cast<double>(validIdx[0]));
                loop_ub = 0;
                i2 = validIdx.size(0);
                b_this.set_size(validIdx.size(0));
                c_this.set_size(validIdx.size(0));
                for (k = 0; k < i2; k++) {
                  if (validIdx[k] != 0) {
                    loop_ub++;
                  }

                  b_this[k] = BoardCoords[((validIdx[k] + BoardCoords.size(0) *
                    j) + BoardCoords.size(0) * BoardCoords.size(1) * (
                    static_cast<int>(coordsToUse[0]) - 1)) - 1];
                  c_this[k] = BoardCoords[((validIdx[k] + BoardCoords.size(0) *
                    j) + BoardCoords.size(0) * BoardCoords.size(1) * (
                    static_cast<int>(coordsToUse[1]) - 1)) - 1];
                }

                if (loop_ub > 5) {
                  i2 = 4;
                } else {
                  i2 = 2;
                }

                polyfit(b_this, c_this, static_cast<double>(i2), currCurve_data,
                        currCurve_size);
                currRad = coordDist / 4.0;
                refCoordValue = BoardCoords[((validIdx[0] + BoardCoords.size(0) *
                  j) + BoardCoords.size(0) * BoardCoords.size(1) * (static_cast<
                  int>(coordsToUse[0]) - 1)) - 1];
                currCoord = currRad + refCoordValue;
                exitg1 = false;
                while ((!exitg1) && (std::abs(currCoord - refCoordValue) <
                                     static_cast<double>(validIdx[0]) * 1.5 *
                                     std::abs(coordDist))) {
                  double currPt[2];
                  boolean_T exitg2;
                  boolean_T p;
                  p = true;
                  k = 0;
                  exitg2 = false;
                  while ((!exitg2) && (k < 2)) {
                    if (!(coordsToUse[k] == static_cast<double>(k) + 1.0)) {
                      p = false;
                      exitg2 = true;
                    } else {
                      k++;
                    }
                  }

                  if (p) {
                    double y;
                    y = currCurve_data[0];
                    i2 = currCurve_size[1];
                    for (k = 0; k <= i2 - 2; k++) {
                      y = currCoord * y + currCurve_data[k + 1];
                    }

                    currPt[0] = currCoord;
                    currPt[1] = y;
                  } else {
                    double y;
                    y = currCurve_data[0];
                    i2 = currCurve_size[1];
                    for (k = 0; k <= i2 - 2; k++) {
                      y = currCoord * y + currCurve_data[k + 1];
                    }

                    currPt[0] = y;
                    currPt[1] = currCoord;
                  }

                  findClosestOnCurve(currPt, std::abs(currRad), currCurve_data,
                                     currCurve_size, coordsToUse, removedIdx,
                                     b_index);
                  if (b_index.size(1) != 0) {
                    newIndices[j] = b_index[0];
                    i2 = removedIdx.size(1);
                    loop_ub = b_index.size(1);
                    removedIdx.set_size(removedIdx.size(0), removedIdx.size(1) +
                                        b_index.size(1));
                    for (k = 0; k < loop_ub; k++) {
                      removedIdx[i2 + k] = b_index[k];
                    }

                    exitg1 = true;
                  } else {
                    currCoord += currRad;
                  }
                }
              }
            }
          }

          //
          // Arguments    : void
          // Return Type  : void
          //
          void Checkerboard::undoLastExpansion()
          {
            Energy = PreviousEnergy;
            switch (static_cast<int>(LastExpandDirection)) {
             case 1:
              {
                int b_this;
                int c_this;
                int i;
                int i1;
                int loop_ub_tmp;
                i = BoardIdx.size(0);
                if (i < 2) {
                  i1 = -1;
                  i = -1;
                } else {
                  i1 = 0;
                  i--;
                }

                b_this = BoardIdx.size(1);
                for (int i2{0}; i2 < b_this; i2++) {
                  loop_ub_tmp = i - i1;
                  for (int i3{0}; i3 < loop_ub_tmp; i3++) {
                    BoardIdx[i3 + loop_ub_tmp * i2] = BoardIdx[((i1 + i3) +
                      BoardIdx.size(0) * i2) + 1];
                  }
                }

                BoardIdx.set_size(i - i1, b_this);
                i = BoardCoords.size(0);
                if (i < 2) {
                  i1 = -1;
                  i = -1;
                } else {
                  i1 = 0;
                  i--;
                }

                b_this = BoardCoords.size(1);
                c_this = BoardCoords.size(2);
                for (int i2{0}; i2 < c_this; i2++) {
                  for (int i3{0}; i3 < b_this; i3++) {
                    loop_ub_tmp = i - i1;
                    for (int i4{0}; i4 < loop_ub_tmp; i4++) {
                      BoardCoords[(i4 + loop_ub_tmp * i3) + loop_ub_tmp * b_this
                        * i2] = BoardCoords[(((i1 + i4) + BoardCoords.size(0) *
                        i3) + BoardCoords.size(0) * BoardCoords.size(1) * i2) +
                        1];
                    }
                  }
                }

                BoardCoords.set_size(i - i1, b_this, c_this);
              }
              break;

             case 2:
              {
                int b_this;
                int c_this;
                int i;
                int loop_ub_tmp;
                i = BoardIdx.size(0) - 2;
                if (i + 1 < 1) {
                  loop_ub_tmp = 0;
                } else {
                  loop_ub_tmp = i + 1;
                }

                b_this = BoardIdx.size(1);
                for (i = 0; i < b_this; i++) {
                  for (int i1{0}; i1 < loop_ub_tmp; i1++) {
                    BoardIdx[i1 + loop_ub_tmp * i] = BoardIdx[i1 + BoardIdx.size
                      (0) * i];
                  }
                }

                BoardIdx.set_size(loop_ub_tmp, b_this);
                i = BoardCoords.size(0) - 2;
                if (i + 1 < 1) {
                  loop_ub_tmp = 0;
                } else {
                  loop_ub_tmp = i + 1;
                }

                b_this = BoardCoords.size(1);
                c_this = BoardCoords.size(2);
                for (i = 0; i < c_this; i++) {
                  for (int i1{0}; i1 < b_this; i1++) {
                    for (int i2{0}; i2 < loop_ub_tmp; i2++) {
                      BoardCoords[(i2 + loop_ub_tmp * i1) + loop_ub_tmp * b_this
                        * i] = BoardCoords[(i2 + BoardCoords.size(0) * i1) +
                        BoardCoords.size(0) * BoardCoords.size(1) * i];
                    }
                  }
                }

                BoardCoords.set_size(loop_ub_tmp, b_this, c_this);
              }
              break;

             case 3:
              {
                int b_this;
                int c_this;
                int i;
                int i1;
                int loop_ub_tmp;
                i = BoardIdx.size(1);
                if (i < 2) {
                  i1 = 0;
                  i = 0;
                } else {
                  i1 = 1;
                }

                b_this = BoardIdx.size(0);
                loop_ub_tmp = i - i1;
                for (i = 0; i < loop_ub_tmp; i++) {
                  for (int i2{0}; i2 < b_this; i2++) {
                    BoardIdx[i2 + b_this * i] = BoardIdx[i2 + BoardIdx.size(0) *
                      (i1 + i)];
                  }
                }

                BoardIdx.set_size(b_this, loop_ub_tmp);
                i = BoardCoords.size(1);
                if (i < 2) {
                  i1 = -1;
                  i = -1;
                } else {
                  i1 = 0;
                  i--;
                }

                b_this = BoardCoords.size(0);
                c_this = BoardCoords.size(2);
                for (int i2{0}; i2 < c_this; i2++) {
                  loop_ub_tmp = i - i1;
                  for (int i3{0}; i3 < loop_ub_tmp; i3++) {
                    for (int i4{0}; i4 < b_this; i4++) {
                      BoardCoords[(i4 + b_this * i3) + b_this * loop_ub_tmp * i2]
                        = BoardCoords[(i4 + BoardCoords.size(0) * ((i1 + i3) + 1))
                        + BoardCoords.size(0) * BoardCoords.size(1) * i2];
                    }
                  }
                }

                BoardCoords.set_size(b_this, i - i1, c_this);
              }
              break;

             case 4:
              {
                int b_this;
                int c_this;
                int i;
                int loop_ub_tmp;
                i = BoardIdx.size(1) - 2;
                if (i + 1 < 1) {
                  loop_ub_tmp = -1;
                } else {
                  loop_ub_tmp = i;
                }

                b_this = BoardIdx.size(0);
                for (i = 0; i <= loop_ub_tmp; i++) {
                  for (int i1{0}; i1 < b_this; i1++) {
                    BoardIdx[i1 + b_this * i] = BoardIdx[i1 + BoardIdx.size(0) *
                      i];
                  }
                }

                BoardIdx.set_size(b_this, loop_ub_tmp + 1);
                i = BoardCoords.size(1) - 2;
                if (i + 1 < 1) {
                  loop_ub_tmp = 0;
                } else {
                  loop_ub_tmp = i + 1;
                }

                b_this = BoardCoords.size(0);
                c_this = BoardCoords.size(2);
                for (i = 0; i < c_this; i++) {
                  for (int i1{0}; i1 < loop_ub_tmp; i1++) {
                    for (int i2{0}; i2 < b_this; i2++) {
                      BoardCoords[(i2 + b_this * i1) + b_this * loop_ub_tmp * i]
                        = BoardCoords[(i2 + BoardCoords.size(0) * i1) +
                        BoardCoords.size(0) * BoardCoords.size(1) * i];
                    }
                  }
                }

                BoardCoords.set_size(b_this, loop_ub_tmp, c_this);
              }
              break;
            }
          }

          //
          // Arguments    : coder::vision::internal::calibration::checkerboard::Checkerboard *in1
          //                const coder::array<float, 2U> &in2
          //                const coder::array<float, 1U> &in3
          //                const coder::array<float, 2U> &in4
          //                const coder::array<float, 2U> &in5
          // Return Type  : void
          //
        }
      }
    }
  }
}

static void b_binary_expand_op(coder::vision::internal::calibration::
  checkerboard::Checkerboard *in1, const coder::array<float, 2U> &in2, const
  coder::array<float, 1U> &in3, const coder::array<float, 2U> &in4, const coder::
  array<float, 2U> &in5)
{
  coder::array<float, 2U> b_in4;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  if (in5.size(1) == 1) {
    stride_0_1 = in4.size(1);
  } else {
    stride_0_1 = in5.size(1);
  }

  b_in4.set_size(1, stride_0_1);
  stride_0_1 = (in4.size(1) != 1);
  stride_1_1 = (in5.size(1) != 1);
  if (in5.size(1) == 1) {
    loop_ub = in4.size(1);
  } else {
    loop_ub = in5.size(1);
  }

  if (static_cast<int>(loop_ub < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      b_in4[i] = in4[i * stride_0_1] + in5[i * stride_1_1];
    }
  } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      b_in4[i] = in4[i * stride_0_1] + in5[i * stride_1_1];
    }
  }

  in1->BoardIdx[in1->BoardIdx.size(0) * 2 + 2] = in1->findNeighbor(in2, in3,
    b_in4);
}

//
// Arguments    : coder::array<float, 2U> &in1
//                const coder::array<float, 2U> &in2
// Return Type  : void
//
static void b_minus(coder::array<float, 2U> &in1, const coder::array<float, 2U>
                    &in2)
{
  coder::array<float, 2U> b_in1;
  int aux_0_1;
  int aux_1_1;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  if (in2.size(1) == 1) {
    i = in1.size(1);
  } else {
    i = in2.size(1);
  }

  b_in1.set_size(3, i);
  stride_0_1 = (in1.size(1) != 1);
  stride_1_1 = (in2.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (in2.size(1) == 1) {
    loop_ub = in1.size(1);
  } else {
    loop_ub = in2.size(1);
  }

  for (i = 0; i < loop_ub; i++) {
    b_in1[3 * i] = in1[3 * aux_0_1] - in2[3 * aux_1_1];
    b_in1[3 * i + 1] = in1[3 * aux_0_1 + 1] - in2[3 * aux_1_1 + 1];
    b_in1[3 * i + 2] = in1[3 * aux_0_1 + 2] - in2[3 * aux_1_1 + 2];
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }

  in1.set_size(3, b_in1.size(1));
  loop_ub = b_in1.size(1);
  if (static_cast<int>(b_in1.size(1) * 3 < 3200)) {
    for (int i1{0}; i1 < loop_ub; i1++) {
      in1[3 * i1] = b_in1[3 * i1];
      in1[3 * i1 + 1] = b_in1[3 * i1 + 1];
      in1[3 * i1 + 2] = b_in1[3 * i1 + 2];
    }
  } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i1 = 0; i1 < loop_ub; i1++) {
      in1[3 * i1] = b_in1[3 * i1];
      in1[3 * i1 + 1] = b_in1[3 * i1 + 1];
      in1[3 * i1 + 2] = b_in1[3 * i1 + 2];
    }
  }
}

//
// Arguments    : coder::vision::internal::calibration::checkerboard::Checkerboard *in1
//                const coder::array<float, 2U> &in2
//                const coder::array<float, 1U> &in3
//                const coder::array<float, 2U> &in4
//                const coder::array<float, 2U> &in5
// Return Type  : void
//
static void binary_expand_op(coder::vision::internal::calibration::checkerboard::
  Checkerboard *in1, const coder::array<float, 2U> &in2, const coder::array<
  float, 1U> &in3, const coder::array<float, 2U> &in4, const coder::array<float,
  2U> &in5)
{
  coder::array<float, 2U> b_in4;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  if (in5.size(1) == 1) {
    stride_0_1 = in4.size(1);
  } else {
    stride_0_1 = in5.size(1);
  }

  b_in4.set_size(1, stride_0_1);
  stride_0_1 = (in4.size(1) != 1);
  stride_1_1 = (in5.size(1) != 1);
  if (in5.size(1) == 1) {
    loop_ub = in4.size(1);
  } else {
    loop_ub = in5.size(1);
  }

  if (static_cast<int>(loop_ub < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      b_in4[i] = in4[i * stride_0_1] + in5[i * stride_1_1];
    }
  } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      b_in4[i] = in4[i * stride_0_1] + in5[i * stride_1_1];
    }
  }

  in1->BoardIdx[in1->BoardIdx.size(0) * 2] = in1->findNeighbor(in2, in3, b_in4);
}

//
// Arguments    : coder::array<double, 2U> &in1
//                const coder::array<double, 2U> &in2
//                const coder::array<double, 2U> &in3
// Return Type  : void
//
static void binary_expand_op(coder::array<double, 2U> &in1, const coder::array<
  double, 2U> &in2, const coder::array<double, 2U> &in3)
{
  int aux_0_1;
  int aux_1_1;
  int i;
  int i1;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  if (in3.size(0) == 1) {
    i = in2.size(0);
  } else {
    i = in3.size(0);
  }

  if (in3.size(1) == 1) {
    i1 = in2.size(1);
  } else {
    i1 = in3.size(1);
  }

  in1.set_size(i, i1);
  stride_0_0 = (in2.size(0) != 1);
  stride_0_1 = (in2.size(1) != 1);
  stride_1_0 = (in3.size(0) != 1);
  stride_1_1 = (in3.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (in3.size(1) == 1) {
    loop_ub = in2.size(1);
  } else {
    loop_ub = in3.size(1);
  }

  for (i = 0; i < loop_ub; i++) {
    int b_loop_ub;
    i1 = in3.size(0);
    if (i1 == 1) {
      b_loop_ub = in2.size(0);
    } else {
      b_loop_ub = i1;
    }

    for (i1 = 0; i1 < b_loop_ub; i1++) {
      in1[i1 + in1.size(0) * i] = (in2[i1 * stride_0_0 + in2.size(0) * aux_0_1]
        + in2[i1 * stride_0_0 + in2.size(0) * aux_0_1]) - in3[i1 * stride_1_0 +
        in3.size(0) * aux_1_1];
    }

    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
}

//
// Arguments    : coder::array<double, 2U> &in1
//                const coder::vision::internal::calibration::checkerboard::Checkerboard *in2
//                const coder::array<int, 1U> &in3
//                const coder::array<int, 1U> &in4
//                const coder::array<int, 1U> &in5
// Return Type  : void
//
static void binary_expand_op(coder::array<double, 2U> &in1, const coder::vision::
  internal::calibration::checkerboard::Checkerboard *in2, const coder::array<int,
  1U> &in3, const coder::array<int, 1U> &in4, const coder::array<int, 1U> &in5)
{
  coder::array<double, 3U> b_in2;
  int b_loop_ub;
  int i1;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  if (in5.size(0) == 1) {
    stride_0_0 = in3.size(0);
  } else {
    stride_0_0 = in5.size(0);
  }

  b_in2.set_size(stride_0_0, 1, in2->BoardCoords.size(2));
  stride_0_0 = (in3.size(0) != 1);
  stride_1_0 = (in5.size(0) != 1);
  loop_ub = in2->BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i1,b_loop_ub)

  for (int i = 0; i < loop_ub; i++) {
    i1 = in5.size(0);
    if (i1 == 1) {
      b_loop_ub = in3.size(0);
    } else {
      b_loop_ub = i1;
    }

    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_in2[i1 + b_in2.size(0) * i] = (in2->BoardCoords[(in3[i1 * stride_0_0] +
        in2->BoardCoords.size(0) * in2->BoardCoords.size(1) * i) - 1] +
        in2->BoardCoords[((in4[i1 * stride_0_0] + in2->BoardCoords.size(0) * 2)
                          + in2->BoardCoords.size(0) * in2->BoardCoords.size(1) *
                          i) - 1]) - 2.0 * in2->BoardCoords[((in5[i1 *
        stride_1_0] + in2->BoardCoords.size(0)) + in2->BoardCoords.size(0) *
        in2->BoardCoords.size(1) * i) - 1];
    }
  }

  coder::b_squeeze(b_in2, in1);
}

//
// Arguments    : coder::array<double, 2U> &in1
//                const coder::vision::internal::calibration::checkerboard::Checkerboard *in2
//                const coder::array<double, 2U> &in3
//                const coder::array<double, 2U> &in4
// Return Type  : void
//
static void binary_expand_op(coder::array<double, 2U> &in1, const coder::vision::
  internal::calibration::checkerboard::Checkerboard *in2, const coder::array<
  double, 2U> &in3, const coder::array<double, 2U> &in4)
{
  coder::array<double, 2U> b_in3;
  int aux_0_1;
  int aux_1_1;
  int i;
  int i1;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  if (in4.size(0) == 1) {
    i = in3.size(0);
  } else {
    i = in4.size(0);
  }

  if (in4.size(1) == 1) {
    i1 = in3.size(1);
  } else {
    i1 = in4.size(1);
  }

  b_in3.set_size(i, i1);
  stride_0_0 = (in3.size(0) != 1);
  stride_0_1 = (in3.size(1) != 1);
  stride_1_0 = (in4.size(0) != 1);
  stride_1_1 = (in4.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (in4.size(1) == 1) {
    loop_ub = in3.size(1);
  } else {
    loop_ub = in4.size(1);
  }

  for (i = 0; i < loop_ub; i++) {
    int b_loop_ub;
    i1 = in4.size(0);
    if (i1 == 1) {
      b_loop_ub = in3.size(0);
    } else {
      b_loop_ub = i1;
    }

    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_in3[i1 + b_in3.size(0) * i] = (in3[i1 * stride_0_0 + in3.size(0) *
        aux_0_1] + in3[i1 * stride_0_0 + in3.size(0) * aux_0_1]) - in4[i1 *
        stride_1_0 + in4.size(0) * aux_1_1];
    }

    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }

  in2->findClosestIndices(b_in3, in1);
}

//
// Arguments    : coder::array<double, 2U> &in1
//                const coder::array<double, 3U> &in2
//                const coder::array<double, 3U> &in3
//                const coder::array<double, 3U> &in4
// Return Type  : void
//
static void binary_expand_op(coder::array<double, 2U> &in1, const coder::array<
  double, 3U> &in2, const coder::array<double, 3U> &in3, const coder::array<
  double, 3U> &in4)
{
  coder::array<double, 3U> b_in2;
  int aux_0_2;
  int aux_1_2;
  int i;
  int i1;
  int loop_ub;
  int stride_0_1;
  int stride_0_2;
  int stride_1_1;
  int stride_1_2;
  if (in4.size(1) == 1) {
    i = in2.size(1);
  } else {
    i = in4.size(1);
  }

  if (in4.size(2) == 1) {
    i1 = in2.size(2);
  } else {
    i1 = in4.size(2);
  }

  b_in2.set_size(1, i, i1);
  stride_0_1 = (in2.size(1) != 1);
  stride_0_2 = (in2.size(2) != 1);
  stride_1_1 = (in4.size(1) != 1);
  stride_1_2 = (in4.size(2) != 1);
  aux_0_2 = 0;
  aux_1_2 = 0;
  if (in4.size(2) == 1) {
    loop_ub = in2.size(2);
  } else {
    loop_ub = in4.size(2);
  }

  for (i = 0; i < loop_ub; i++) {
    int b_loop_ub;
    i1 = in4.size(1);
    if (i1 == 1) {
      b_loop_ub = in2.size(1);
    } else {
      b_loop_ub = i1;
    }

    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_in2[i1 + b_in2.size(1) * i] = (in2[i1 * stride_0_1 + in2.size(1) *
        aux_0_2] + in3[i1 * stride_0_1 + in3.size(1) * aux_0_2]) - 2.0 * in4[i1 *
        stride_1_1 + in4.size(1) * aux_1_2];
    }

    aux_1_2 += stride_1_2;
    aux_0_2 += stride_0_2;
  }

  coder::squeeze(b_in2, in1);
}

//
// Arguments    : coder::array<float, 1U> &in1
//                const coder::array<float, 1U> &in2
//                const coder::array<float, 1U> &in3
// Return Type  : void
//
static void binary_expand_op(coder::array<float, 1U> &in1, const coder::array<
  float, 1U> &in2, const coder::array<float, 1U> &in3)
{
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  if (in3.size(0) == 1) {
    loop_ub = in2.size(0);
  } else {
    loop_ub = in3.size(0);
  }

  if (loop_ub == 1) {
    loop_ub = in2.size(0);
  } else if (in3.size(0) == 1) {
    loop_ub = in2.size(0);
  } else {
    loop_ub = in3.size(0);
  }

  in1.set_size(loop_ub);
  stride_0_0 = (in2.size(0) != 1);
  stride_1_0 = (in2.size(0) != 1);
  stride_2_0 = (in3.size(0) != 1);
  if (in3.size(0) == 1) {
    loop_ub = in2.size(0);
  } else {
    loop_ub = in3.size(0);
  }

  if (loop_ub == 1) {
    loop_ub = in2.size(0);
  } else if (in3.size(0) == 1) {
    loop_ub = in2.size(0);
  } else {
    loop_ub = in3.size(0);
  }

  if (static_cast<int>(loop_ub < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      in1[i] = in2[i * stride_0_0] + 1.5F * in2[i * stride_1_0] * (1.0F - in3[i *
        stride_2_0]);
    }
  } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      in1[i] = in2[i * stride_0_0] + 1.5F * in2[i * stride_1_0] * (1.0F - in3[i *
        stride_2_0]);
    }
  }
}

//
// Arguments    : coder::array<double, 2U> &in1
//                const coder::vision::internal::calibration::checkerboard::Checkerboard *in2
//                const coder::array<int, 1U> &in3
//                const coder::array<double, 2U> &in4
//                const coder::array<int, 1U> &in5
//                const coder::array<int, 1U> &in6
// Return Type  : void
//
static void binary_expand_op(coder::array<double, 2U> &in1, const coder::vision::
  internal::calibration::checkerboard::Checkerboard *in2, const coder::array<int,
  1U> &in3, const coder::array<double, 2U> &in4, const coder::array<int, 1U>
  &in5, const coder::array<int, 1U> &in6)
{
  coder::array<double, 3U> b_in2;
  int b_in4;
  int b_loop_ub;
  int c_in4;
  int d_in4;
  int i1;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  b_in4 = static_cast<int>(in4[0]);
  c_in4 = static_cast<int>(in4[2]);
  d_in4 = static_cast<int>(in4[1]);
  if (in6.size(0) == 1) {
    stride_0_0 = in3.size(0);
  } else {
    stride_0_0 = in6.size(0);
  }

  b_in2.set_size(stride_0_0, 1, in2->BoardCoords.size(2));
  stride_0_0 = (in3.size(0) != 1);
  stride_1_0 = (in6.size(0) != 1);
  loop_ub = in2->BoardCoords.size(2);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(i1,b_loop_ub)

  for (int i = 0; i < loop_ub; i++) {
    i1 = in6.size(0);
    if (i1 == 1) {
      b_loop_ub = in3.size(0);
    } else {
      b_loop_ub = i1;
    }

    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_in2[i1 + b_in2.size(0) * i] = (in2->BoardCoords[((in3[i1 * stride_0_0] +
        in2->BoardCoords.size(0) * (b_in4 - 1)) + in2->BoardCoords.size(0) *
        in2->BoardCoords.size(1) * i) - 1] + in2->BoardCoords[((in5[i1 *
        stride_0_0] + in2->BoardCoords.size(0) * (c_in4 - 1)) +
        in2->BoardCoords.size(0) * in2->BoardCoords.size(1) * i) - 1]) - 2.0 *
        in2->BoardCoords[((in6[i1 * stride_1_0] + in2->BoardCoords.size(0) *
                           (d_in4 - 1)) + in2->BoardCoords.size(0) *
                          in2->BoardCoords.size(1) * i) - 1];
    }
  }

  coder::b_squeeze(b_in2, in1);
}

//
// Arguments    : coder::vision::internal::calibration::checkerboard::Checkerboard *in1
//                const coder::array<float, 2U> &in2
//                const coder::array<float, 1U> &in3
//                const coder::array<float, 2U> &in4
//                const coder::array<float, 2U> &in5
// Return Type  : void
//
static void c_binary_expand_op(coder::vision::internal::calibration::
  checkerboard::Checkerboard *in1, const coder::array<float, 2U> &in2, const
  coder::array<float, 1U> &in3, const coder::array<float, 2U> &in4, const coder::
  array<float, 2U> &in5)
{
  coder::array<float, 2U> b_in4;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  if (in5.size(1) == 1) {
    stride_0_1 = in4.size(1);
  } else {
    stride_0_1 = in5.size(1);
  }

  b_in4.set_size(1, stride_0_1);
  stride_0_1 = (in4.size(1) != 1);
  stride_1_1 = (in5.size(1) != 1);
  if (in5.size(1) == 1) {
    loop_ub = in4.size(1);
  } else {
    loop_ub = in5.size(1);
  }

  if (static_cast<int>(loop_ub < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      b_in4[i] = in4[i * stride_0_1] + in5[i * stride_1_1];
    }
  } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      b_in4[i] = in4[i * stride_0_1] + in5[i * stride_1_1];
    }
  }

  in1->BoardIdx[2] = in1->findNeighbor(in2, in3, b_in4);
}

//
// Arguments    : coder::array<float, 2U> &in1
//                const coder::array<float, 2U> &in2
//                const coder::array<float, 2U> &in3
// Return Type  : void
//
static void c_binary_expand_op(coder::array<float, 2U> &in1, const coder::array<
  float, 2U> &in2, const coder::array<float, 2U> &in3)
{
  coder::array<float, 2U> b_in2;
  int aux_0_1;
  int aux_1_1;
  int aux_2_1;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  int stride_2_1;
  if (in1.size(1) == 1) {
    if (in3.size(1) == 1) {
      i = in2.size(1);
    } else {
      i = in3.size(1);
    }
  } else {
    i = in1.size(1);
  }

  b_in2.set_size(3, i);
  stride_0_1 = (in2.size(1) != 1);
  stride_1_1 = (in3.size(1) != 1);
  stride_2_1 = (in1.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  aux_2_1 = 0;
  if (in1.size(1) == 1) {
    if (in3.size(1) == 1) {
      loop_ub = in2.size(1);
    } else {
      loop_ub = in3.size(1);
    }
  } else {
    loop_ub = in1.size(1);
  }

  for (i = 0; i < loop_ub; i++) {
    b_in2[3 * i] = (in2[3 * aux_0_1] + in3[3 * aux_1_1]) - 2.0F * in1[3 *
      aux_2_1];
    b_in2[3 * i + 1] = (in2[3 * aux_0_1 + 1] + in3[3 * aux_1_1 + 1]) - 2.0F *
      in1[3 * aux_2_1 + 1];
    b_in2[3 * i + 2] = (in2[3 * aux_0_1 + 2] + in3[3 * aux_1_1 + 2]) - 2.0F *
      in1[3 * aux_2_1 + 2];
    aux_2_1 += stride_2_1;
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }

  in1.set_size(3, b_in2.size(1));
  loop_ub = b_in2.size(1);
  if (static_cast<int>(b_in2.size(1) * 3 < 3200)) {
    for (int i1{0}; i1 < loop_ub; i1++) {
      in1[3 * i1] = b_in2[3 * i1];
      in1[3 * i1 + 1] = b_in2[3 * i1 + 1];
      in1[3 * i1 + 2] = b_in2[3 * i1 + 2];
    }
  } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i1 = 0; i1 < loop_ub; i1++) {
      in1[3 * i1] = b_in2[3 * i1];
      in1[3 * i1 + 1] = b_in2[3 * i1 + 1];
      in1[3 * i1 + 2] = b_in2[3 * i1 + 2];
    }
  }
}

//
// Arguments    : coder::vision::internal::calibration::checkerboard::Checkerboard *in1
//                const coder::array<float, 2U> &in2
//                const coder::array<float, 1U> &in3
//                const coder::array<float, 2U> &in4
//                const coder::array<float, 2U> &in5
// Return Type  : void
//
static void d_binary_expand_op(coder::vision::internal::calibration::
  checkerboard::Checkerboard *in1, const coder::array<float, 2U> &in2, const
  coder::array<float, 1U> &in3, const coder::array<float, 2U> &in4, const coder::
  array<float, 2U> &in5)
{
  coder::array<float, 2U> b_in4;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  if (in5.size(1) == 1) {
    stride_0_1 = in4.size(1);
  } else {
    stride_0_1 = in5.size(1);
  }

  b_in4.set_size(1, stride_0_1);
  stride_0_1 = (in4.size(1) != 1);
  stride_1_1 = (in5.size(1) != 1);
  if (in5.size(1) == 1) {
    loop_ub = in4.size(1);
  } else {
    loop_ub = in5.size(1);
  }

  if (static_cast<int>(loop_ub < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      b_in4[i] = in4[i * stride_0_1] + in5[i * stride_1_1];
    }
  } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      b_in4[i] = in4[i * stride_0_1] + in5[i * stride_1_1];
    }
  }

  in1->BoardIdx[0] = in1->findNeighbor(in2, in3, b_in4);
}

//
// Arguments    : coder::array<boolean_T, 1U> &in1
//                const coder::array<double, 2U> &in2
//                const coder::array<double, 2U> &in3
// Return Type  : void
//
static void e_binary_expand_op(coder::array<boolean_T, 1U> &in1, const coder::
  array<double, 2U> &in2, const coder::array<double, 2U> &in3)
{
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  if (in3.size(0) == 1) {
    stride_0_0 = in2.size(0);
  } else {
    stride_0_0 = in3.size(0);
  }

  in1.set_size(stride_0_0);
  stride_0_0 = (in2.size(0) != 1);
  stride_1_0 = (in3.size(0) != 1);
  if (in3.size(0) == 1) {
    loop_ub = in2.size(0);
  } else {
    loop_ub = in3.size(0);
  }

  if (static_cast<int>(loop_ub < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      in1[i] = ((in2[i * stride_0_0] == 0.0) || (in3[i * stride_1_0] == 0.0));
    }
  } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      in1[i] = ((in2[i * stride_0_0] == 0.0) || (in3[i * stride_1_0] == 0.0));
    }
  }
}

//
// Arguments    : coder::array<float, 2U> &in1
//                const coder::array<float, 2U> &in2
// Return Type  : void
//
static void minus(coder::array<float, 2U> &in1, const coder::array<float, 2U>
                  &in2)
{
  coder::array<float, 2U> b_in1;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  if (in2.size(1) == 1) {
    stride_0_1 = in1.size(1);
  } else {
    stride_0_1 = in2.size(1);
  }

  b_in1.set_size(1, stride_0_1);
  stride_0_1 = (in1.size(1) != 1);
  stride_1_1 = (in2.size(1) != 1);
  if (in2.size(1) == 1) {
    loop_ub = in1.size(1);
  } else {
    loop_ub = in2.size(1);
  }

  if (static_cast<int>(loop_ub < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      b_in1[i] = in1[i * stride_0_1] - in2[i * stride_1_1];
    }
  } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      b_in1[i] = in1[i * stride_0_1] - in2[i * stride_1_1];
    }
  }

  in1.set_size(1, b_in1.size(1));
  loop_ub = b_in1.size(1);
  if (static_cast<int>(b_in1.size(1) < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      in1[i] = b_in1[i];
    }
  } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      in1[i] = b_in1[i];
    }
  }
}

//
// Arguments    : void
// Return Type  : boolean_T
//
namespace coder
{
  namespace vision
  {
    namespace internal
    {
      namespace calibration
      {
        namespace checkerboard
        {
          boolean_T Checkerboard::b_expandBoardOnce()
          {
            array<double, 3U> b_this;
            array<double, 3U> c_this;
            array<double, 3U> r;
            array<double, 2U> b_index;
            array<double, 2U> idx;
            array<double, 2U> newIndices;
            array<double, 2U> p1;
            array<double, 2U> p2;
            array<double, 2U> predictedPoints;
            array<double, 2U> removedIdx;
            array<double, 1U> d_this;
            array<double, 1U> e_this;
            array<double, 1U> validIdx;
            array<int, 1U> r1;
            array<int, 1U> r2;
            array<int, 1U> r3;
            array<int, 1U> r4;
            array<int, 1U> r5;
            array<boolean_T, 1U> invalidPts;
            double currCurve_data[5];
            double coordDist;
            double moveDistMultiplier;
            double refCoordValue;
            int i;
            boolean_T success;
            PreviousEnergy = Energy;
            i = 0;
            int exitg1;
            do {
              exitg1 = 0;
              if (i < 4) {
                if (!IsDirectionBad[i]) {
                  float oldEnergy;
                  LastExpandDirection = static_cast<double>(i) + 1.0;
                  oldEnergy = (Energy + static_cast<float>(BoardIdx.size(0) *
                    BoardIdx.size(1))) / static_cast<float>(BoardIdx.size(0) *
                    BoardIdx.size(1));
                  switch (i + 1) {
                   case 1:
                    {
                      int loop_ub;
                      if (IsDistortionHigh) {
                        c_fitPolynomialIndices(newIndices);
                      } else {
                        int b_loop_ub;
                        int numCols;
                        b_this.set_size(1, BoardCoords.size(1), BoardCoords.size
                                        (2));
                        loop_ub = BoardCoords.size(2);
                        for (int b_i{0}; b_i < loop_ub; b_i++) {
                          b_loop_ub = BoardCoords.size(1);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            b_this[i1 + b_this.size(1) * b_i] = BoardCoords
                              [(BoardCoords.size(0) * i1 + BoardCoords.size(0) *
                                BoardCoords.size(1) * b_i) + 1];
                          }
                        }

                        squeeze(b_this, p1);
                        b_this.set_size(1, BoardCoords.size(1), BoardCoords.size
                                        (2));
                        loop_ub = BoardCoords.size(2);
                        for (int b_i{0}; b_i < loop_ub; b_i++) {
                          b_loop_ub = BoardCoords.size(1);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            b_this[i1 + b_this.size(1) * b_i] =
                              BoardCoords[BoardCoords.size(0) * i1 +
                              BoardCoords.size(0) * BoardCoords.size(1) * b_i];
                          }
                        }

                        squeeze(b_this, p2);
                        if ((p2.size(0) == p1.size(0)) && (p2.size(1) == p1.size
                             (1))) {
                          predictedPoints.set_size(p2.size(0), p2.size(1));
                          loop_ub = p2.size(0) * p2.size(1);
                          for (int b_i{0}; b_i < loop_ub; b_i++) {
                            predictedPoints[b_i] = (p2[b_i] + p2[b_i]) - p1[b_i];
                          }
                        } else {
                          binary_expand_op(predictedPoints, p2, p1);
                        }

                        if (p1.size(0) == p2.size(0)) {
                          invalidPts.set_size(p1.size(0));
                          loop_ub = p1.size(0);
                          for (int b_i{0}; b_i < loop_ub; b_i++) {
                            invalidPts[b_i] = ((p1[b_i] == 0.0) || (p2[b_i] ==
                              0.0));
                          }
                        } else {
                          e_binary_expand_op(invalidPts, p1, p2);
                        }

                        b_loop_ub = invalidPts.size(0) - 1;
                        numCols = 0;
                        for (loop_ub = 0; loop_ub <= b_loop_ub; loop_ub++) {
                          if (invalidPts[loop_ub]) {
                            numCols++;
                          }
                        }

                        r2.set_size(numCols);
                        numCols = 0;
                        for (loop_ub = 0; loop_ub <= b_loop_ub; loop_ub++) {
                          if (invalidPts[loop_ub]) {
                            r2[numCols] = loop_ub + 1;
                            numCols++;
                          }
                        }

                        loop_ub = predictedPoints.size(1);
                        for (int b_i{0}; b_i < loop_ub; b_i++) {
                          b_loop_ub = r2.size(0);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            predictedPoints[(r2[i1] + predictedPoints.size(0) *
                                             b_i) - 1] = rtNaN;
                          }
                        }

                        b_findClosestIndices(predictedPoints, newIndices);
                      }

                      expandBoardUp(newIndices, p1, r);
                      BoardIdx.set_size(p1.size(0), p1.size(1));
                      loop_ub = p1.size(0) * p1.size(1);
                      for (int b_i{0}; b_i < loop_ub; b_i++) {
                        BoardIdx[b_i] = p1[b_i];
                      }

                      BoardCoords.set_size(r.size(0), r.size(1), r.size(2));
                      loop_ub = r.size(0) * r.size(1) * r.size(2);
                      for (int b_i{0}; b_i < loop_ub; b_i++) {
                        BoardCoords[b_i] = r[b_i];
                      }

                      oldEnergy = computeNewEnergyVertical(oldEnergy);
                    }
                    break;

                   case 2:
                    {
                      int b_i;
                      int loop_ub;
                      int numCols;
                      numCols = BoardCoords.size(0);
                      if (numCols < numCols - 2) {
                        idx.set_size(1, 0);
                      } else {
                        idx.set_size(1, 3);
                        for (b_i = 0; b_i < 3; b_i++) {
                          idx[b_i] = numCols - b_i;
                        }
                      }

                      if (IsDistortionHigh) {
                        double coordsToUse[2];
                        int b_loop_ub;
                        findIndependentVar(idx, coordsToUse);
                        newIndices.set_size(1, BoardCoords.size(1));
                        loop_ub = BoardCoords.size(1);
                        for (b_i = 0; b_i < loop_ub; b_i++) {
                          newIndices[b_i] = 0.0;
                        }

                        removedIdx.set_size(1, 0);
                        b_i = newIndices.size(1);
                        for (int j{0}; j < b_i; j++) {
                          int i1;
                          validIdx.set_size(BoardCoords.size(0));
                          loop_ub = BoardCoords.size(0);
                          for (i1 = 0; i1 < loop_ub; i1++) {
                            validIdx[i1] = BoardCoords[(i1 + BoardCoords.size(0)
                              * j) + BoardCoords.size(0) * BoardCoords.size(1) *
                              (static_cast<int>(coordsToUse[0]) - 1)];
                          }

                          eml_find(validIdx, r1);
                          validIdx.set_size(r1.size(0));
                          loop_ub = r1.size(0);
                          for (i1 = 0; i1 < loop_ub; i1++) {
                            validIdx[i1] = r1[i1];
                          }

                          if (validIdx.size(0) >= 2) {
                            double currCoord;
                            double currRad;
                            int currCurve_size[2];
                            boolean_T exitg2;
                            findSearchParams(idx, validIdx, static_cast<double>
                                             (j) + 1.0, coordsToUse, &coordDist,
                                             &moveDistMultiplier, &refCoordValue);
                            numCols = 0;
                            i1 = validIdx.size(0);
                            for (b_loop_ub = 0; b_loop_ub < i1; b_loop_ub++) {
                              if (validIdx[b_loop_ub] != 0.0) {
                                numCols++;
                              }
                            }

                            d_this.set_size(validIdx.size(0));
                            loop_ub = validIdx.size(0);
                            for (i1 = 0; i1 < loop_ub; i1++) {
                              d_this[i1] = BoardCoords[((static_cast<int>
                                (validIdx[i1]) + BoardCoords.size(0) * j) +
                                BoardCoords.size(0) * BoardCoords.size(1) * (
                                static_cast<int>(coordsToUse[0]) - 1)) - 1];
                            }

                            e_this.set_size(validIdx.size(0));
                            loop_ub = validIdx.size(0);
                            for (i1 = 0; i1 < loop_ub; i1++) {
                              e_this[i1] = BoardCoords[((static_cast<int>
                                (validIdx[i1]) + BoardCoords.size(0) * j) +
                                BoardCoords.size(0) * BoardCoords.size(1) * (
                                static_cast<int>(coordsToUse[1]) - 1)) - 1];
                            }

                            if (numCols > 5) {
                              i1 = 4;
                            } else {
                              i1 = 2;
                            }

                            polyfit(d_this, e_this, static_cast<double>(i1),
                                    currCurve_data, currCurve_size);
                            currRad = coordDist / 4.0;
                            refCoordValue = BoardCoords[((static_cast<int>
                              (refCoordValue) + BoardCoords.size(0) * j) +
                              BoardCoords.size(0) * BoardCoords.size(1) * (
                              static_cast<int>(coordsToUse[0]) - 1)) - 1];
                            currCoord = currRad + refCoordValue;
                            exitg2 = false;
                            while ((!exitg2) && (std::abs(currCoord -
                                     refCoordValue) < moveDistMultiplier * 1.5 *
                                                 std::abs(coordDist))) {
                              double currPt[2];
                              boolean_T exitg3;
                              boolean_T p;
                              p = true;
                              b_loop_ub = 0;
                              exitg3 = false;
                              while ((!exitg3) && (b_loop_ub < 2)) {
                                if (!(coordsToUse[b_loop_ub] == static_cast<
                                      double>(b_loop_ub) + 1.0)) {
                                  p = false;
                                  exitg3 = true;
                                } else {
                                  b_loop_ub++;
                                }
                              }

                              if (p) {
                                double y;
                                y = currCurve_data[0];
                                i1 = currCurve_size[1];
                                for (b_loop_ub = 0; b_loop_ub <= i1 - 2;
                                     b_loop_ub++) {
                                  y = currCoord * y + currCurve_data[b_loop_ub +
                                    1];
                                }

                                currPt[0] = currCoord;
                                currPt[1] = y;
                              } else {
                                double y;
                                y = currCurve_data[0];
                                i1 = currCurve_size[1];
                                for (b_loop_ub = 0; b_loop_ub <= i1 - 2;
                                     b_loop_ub++) {
                                  y = currCoord * y + currCurve_data[b_loop_ub +
                                    1];
                                }

                                currPt[0] = y;
                                currPt[1] = currCoord;
                              }

                              findClosestOnCurve(currPt, std::abs(currRad),
                                                 currCurve_data, currCurve_size,
                                                 coordsToUse, removedIdx,
                                                 b_index);
                              if (b_index.size(1) != 0) {
                                newIndices[j] = b_index[0];
                                i1 = removedIdx.size(1);
                                loop_ub = b_index.size(1);
                                removedIdx.set_size(removedIdx.size(0),
                                                    removedIdx.size(1) +
                                                    b_index.size(1));
                                for (numCols = 0; numCols < loop_ub; numCols++)
                                {
                                  removedIdx[i1 + numCols] = b_index[numCols];
                                }

                                exitg2 = true;
                              } else {
                                currCoord += currRad;
                              }
                            }
                          }
                        }

                        numCols = 0;
                        b_i = newIndices.size(1);
                        for (b_loop_ub = 0; b_loop_ub < b_i; b_loop_ub++) {
                          if (newIndices[b_loop_ub] != 0.0) {
                            numCols++;
                          }
                        }

                        if (numCols < 4) {
                          b_loop_ub = newIndices.size(1);
                          for (loop_ub = 0; loop_ub < b_loop_ub; loop_ub++) {
                            if (newIndices[loop_ub] > 0.0) {
                              newIndices[loop_ub] = 0.0;
                            }
                          }
                        }
                      } else {
                        int b_loop_ub;
                        numCols = static_cast<int>(idx[1]);
                        b_this.set_size(1, BoardCoords.size(1), BoardCoords.size
                                        (2));
                        loop_ub = BoardCoords.size(2);
                        for (b_i = 0; b_i < loop_ub; b_i++) {
                          b_loop_ub = BoardCoords.size(1);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            b_this[i1 + b_this.size(1) * b_i] = BoardCoords
                              [((numCols + BoardCoords.size(0) * i1) +
                                BoardCoords.size(0) * BoardCoords.size(1) * b_i)
                              - 1];
                          }
                        }

                        squeeze(b_this, p1);
                        numCols = static_cast<int>(idx[0]);
                        b_this.set_size(1, BoardCoords.size(1), BoardCoords.size
                                        (2));
                        loop_ub = BoardCoords.size(2);
                        for (b_i = 0; b_i < loop_ub; b_i++) {
                          b_loop_ub = BoardCoords.size(1);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            b_this[i1 + b_this.size(1) * b_i] = BoardCoords
                              [((numCols + BoardCoords.size(0) * i1) +
                                BoardCoords.size(0) * BoardCoords.size(1) * b_i)
                              - 1];
                          }
                        }

                        squeeze(b_this, p2);
                        if ((p2.size(0) == p1.size(0)) && (p2.size(1) == p1.size
                             (1))) {
                          predictedPoints.set_size(p2.size(0), p2.size(1));
                          loop_ub = p2.size(0) * p2.size(1);
                          for (b_i = 0; b_i < loop_ub; b_i++) {
                            predictedPoints[b_i] = (p2[b_i] + p2[b_i]) - p1[b_i];
                          }
                        } else {
                          binary_expand_op(predictedPoints, p2, p1);
                        }

                        if (p1.size(0) == p2.size(0)) {
                          invalidPts.set_size(p1.size(0));
                          loop_ub = p1.size(0);
                          for (b_i = 0; b_i < loop_ub; b_i++) {
                            invalidPts[b_i] = ((p1[b_i] == 0.0) || (p2[b_i] ==
                              0.0));
                          }
                        } else {
                          e_binary_expand_op(invalidPts, p1, p2);
                        }

                        b_loop_ub = invalidPts.size(0) - 1;
                        numCols = 0;
                        for (loop_ub = 0; loop_ub <= b_loop_ub; loop_ub++) {
                          if (invalidPts[loop_ub]) {
                            numCols++;
                          }
                        }

                        r5.set_size(numCols);
                        numCols = 0;
                        for (loop_ub = 0; loop_ub <= b_loop_ub; loop_ub++) {
                          if (invalidPts[loop_ub]) {
                            r5[numCols] = loop_ub + 1;
                            numCols++;
                          }
                        }

                        loop_ub = predictedPoints.size(1);
                        for (b_i = 0; b_i < loop_ub; b_i++) {
                          b_loop_ub = r5.size(0);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            predictedPoints[(r5[i1] + predictedPoints.size(0) *
                                             b_i) - 1] = rtNaN;
                          }
                        }

                        b_findClosestIndices(predictedPoints, newIndices);
                      }

                      expandBoardDown(newIndices, p1, r);
                      BoardIdx.set_size(p1.size(0), p1.size(1));
                      loop_ub = p1.size(0) * p1.size(1);
                      for (b_i = 0; b_i < loop_ub; b_i++) {
                        BoardIdx[b_i] = p1[b_i];
                      }

                      BoardCoords.set_size(r.size(0), r.size(1), r.size(2));
                      loop_ub = r.size(0) * r.size(1) * r.size(2);
                      for (b_i = 0; b_i < loop_ub; b_i++) {
                        BoardCoords[b_i] = r[b_i];
                      }

                      idx.set_size(1, idx.size(1));
                      loop_ub = idx.size(1) - 1;
                      for (b_i = 0; b_i <= loop_ub; b_i++) {
                        idx[b_i] = idx[b_i] + 1.0;
                      }

                      oldEnergy = computeNewEnergyVertical(idx, oldEnergy);
                    }
                    break;

                   case 3:
                    {
                      int loop_ub;
                      if (IsDistortionHigh) {
                        d_fitPolynomialIndices(newIndices);
                      } else {
                        int b_loop_ub;
                        int numCols;
                        c_this.set_size(BoardCoords.size(0), 1, BoardCoords.size
                                        (2));
                        loop_ub = BoardCoords.size(2);
                        for (int b_i{0}; b_i < loop_ub; b_i++) {
                          b_loop_ub = BoardCoords.size(0);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            c_this[i1 + c_this.size(0) * b_i] = BoardCoords[(i1
                              + BoardCoords.size(0)) + BoardCoords.size(0) *
                              BoardCoords.size(1) * b_i];
                          }
                        }

                        b_squeeze(c_this, p1);
                        c_this.set_size(BoardCoords.size(0), 1, BoardCoords.size
                                        (2));
                        loop_ub = BoardCoords.size(2);
                        for (int b_i{0}; b_i < loop_ub; b_i++) {
                          b_loop_ub = BoardCoords.size(0);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            c_this[i1 + c_this.size(0) * b_i] = BoardCoords[i1 +
                              BoardCoords.size(0) * BoardCoords.size(1) * b_i];
                          }
                        }

                        b_squeeze(c_this, p2);
                        if ((p2.size(0) == p1.size(0)) && (p2.size(1) == p1.size
                             (1))) {
                          predictedPoints.set_size(p2.size(0), p2.size(1));
                          loop_ub = p2.size(0) * p2.size(1);
                          for (int b_i{0}; b_i < loop_ub; b_i++) {
                            predictedPoints[b_i] = (p2[b_i] + p2[b_i]) - p1[b_i];
                          }
                        } else {
                          binary_expand_op(predictedPoints, p2, p1);
                        }

                        if (p1.size(0) == p2.size(0)) {
                          invalidPts.set_size(p1.size(0));
                          loop_ub = p1.size(0);
                          for (int b_i{0}; b_i < loop_ub; b_i++) {
                            invalidPts[b_i] = ((p1[b_i] == 0.0) || (p2[b_i] ==
                              0.0));
                          }
                        } else {
                          e_binary_expand_op(invalidPts, p1, p2);
                        }

                        b_loop_ub = invalidPts.size(0) - 1;
                        numCols = 0;
                        for (loop_ub = 0; loop_ub <= b_loop_ub; loop_ub++) {
                          if (invalidPts[loop_ub]) {
                            numCols++;
                          }
                        }

                        r3.set_size(numCols);
                        numCols = 0;
                        for (loop_ub = 0; loop_ub <= b_loop_ub; loop_ub++) {
                          if (invalidPts[loop_ub]) {
                            r3[numCols] = loop_ub + 1;
                            numCols++;
                          }
                        }

                        loop_ub = predictedPoints.size(1);
                        for (int b_i{0}; b_i < loop_ub; b_i++) {
                          b_loop_ub = r3.size(0);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            predictedPoints[(r3[i1] + predictedPoints.size(0) *
                                             b_i) - 1] = rtNaN;
                          }
                        }

                        b_findClosestIndices(predictedPoints, newIndices);
                      }

                      expandBoardLeft(newIndices, p1, r);
                      BoardIdx.set_size(p1.size(0), p1.size(1));
                      loop_ub = p1.size(0) * p1.size(1);
                      for (int b_i{0}; b_i < loop_ub; b_i++) {
                        BoardIdx[b_i] = p1[b_i];
                      }

                      BoardCoords.set_size(r.size(0), r.size(1), r.size(2));
                      loop_ub = r.size(0) * r.size(1) * r.size(2);
                      for (int b_i{0}; b_i < loop_ub; b_i++) {
                        BoardCoords[b_i] = r[b_i];
                      }

                      oldEnergy = computeNewEnergyHorizontal(oldEnergy);
                    }
                    break;

                   default:
                    {
                      int loop_ub;
                      int numCols;
                      numCols = BoardCoords.size(1);
                      if (numCols < numCols - 2) {
                        idx.set_size(1, 0);
                      } else {
                        idx.set_size(1, 3);
                        for (int b_i{0}; b_i < 3; b_i++) {
                          idx[b_i] = numCols - b_i;
                        }
                      }

                      if (IsDistortionHigh) {
                        b_fitPolynomialIndices(idx, newIndices);
                      } else {
                        int b_loop_ub;
                        numCols = static_cast<int>(idx[1]);
                        c_this.set_size(BoardCoords.size(0), 1, BoardCoords.size
                                        (2));
                        loop_ub = BoardCoords.size(2);
                        for (int b_i{0}; b_i < loop_ub; b_i++) {
                          b_loop_ub = BoardCoords.size(0);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            c_this[i1 + c_this.size(0) * b_i] = BoardCoords[(i1
                              + BoardCoords.size(0) * (numCols - 1)) +
                              BoardCoords.size(0) * BoardCoords.size(1) * b_i];
                          }
                        }

                        b_squeeze(c_this, p1);
                        numCols = static_cast<int>(idx[0]);
                        c_this.set_size(BoardCoords.size(0), 1, BoardCoords.size
                                        (2));
                        loop_ub = BoardCoords.size(2);
                        for (int b_i{0}; b_i < loop_ub; b_i++) {
                          b_loop_ub = BoardCoords.size(0);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            c_this[i1 + c_this.size(0) * b_i] = BoardCoords[(i1
                              + BoardCoords.size(0) * (numCols - 1)) +
                              BoardCoords.size(0) * BoardCoords.size(1) * b_i];
                          }
                        }

                        b_squeeze(c_this, p2);
                        if ((p2.size(0) == p1.size(0)) && (p2.size(1) == p1.size
                             (1))) {
                          predictedPoints.set_size(p2.size(0), p2.size(1));
                          loop_ub = p2.size(0) * p2.size(1);
                          for (int b_i{0}; b_i < loop_ub; b_i++) {
                            predictedPoints[b_i] = (p2[b_i] + p2[b_i]) - p1[b_i];
                          }
                        } else {
                          binary_expand_op(predictedPoints, p2, p1);
                        }

                        if (p1.size(0) == p2.size(0)) {
                          invalidPts.set_size(p1.size(0));
                          loop_ub = p1.size(0);
                          for (int b_i{0}; b_i < loop_ub; b_i++) {
                            invalidPts[b_i] = ((p1[b_i] == 0.0) || (p2[b_i] ==
                              0.0));
                          }
                        } else {
                          e_binary_expand_op(invalidPts, p1, p2);
                        }

                        b_loop_ub = invalidPts.size(0) - 1;
                        numCols = 0;
                        for (loop_ub = 0; loop_ub <= b_loop_ub; loop_ub++) {
                          if (invalidPts[loop_ub]) {
                            numCols++;
                          }
                        }

                        r4.set_size(numCols);
                        numCols = 0;
                        for (loop_ub = 0; loop_ub <= b_loop_ub; loop_ub++) {
                          if (invalidPts[loop_ub]) {
                            r4[numCols] = loop_ub + 1;
                            numCols++;
                          }
                        }

                        loop_ub = predictedPoints.size(1);
                        for (int b_i{0}; b_i < loop_ub; b_i++) {
                          b_loop_ub = r4.size(0);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            predictedPoints[(r4[i1] + predictedPoints.size(0) *
                                             b_i) - 1] = rtNaN;
                          }
                        }

                        b_findClosestIndices(predictedPoints, newIndices);
                      }

                      expandBoardRight(newIndices, p1, r);
                      BoardIdx.set_size(p1.size(0), p1.size(1));
                      loop_ub = p1.size(0) * p1.size(1);
                      for (int b_i{0}; b_i < loop_ub; b_i++) {
                        BoardIdx[b_i] = p1[b_i];
                      }

                      BoardCoords.set_size(r.size(0), r.size(1), r.size(2));
                      loop_ub = r.size(0) * r.size(1) * r.size(2);
                      for (int b_i{0}; b_i < loop_ub; b_i++) {
                        BoardCoords[b_i] = r[b_i];
                      }

                      idx.set_size(1, idx.size(1));
                      loop_ub = idx.size(1) - 1;
                      for (int b_i{0}; b_i <= loop_ub; b_i++) {
                        idx[b_i] = idx[b_i] + 1.0;
                      }

                      oldEnergy = computeNewEnergyHorizontal(idx, oldEnergy);
                    }
                    break;
                  }

                  Energy = oldEnergy;
                  if (Energy < PreviousEnergy) {
                    success = true;
                    exitg1 = 1;
                  } else {
                    undoLastExpansion();
                    IsDirectionBad[i] = true;
                    i++;
                  }
                } else {
                  i++;
                }
              } else {
                success = false;
                exitg1 = 1;
              }
            } while (exitg1 == 0);

            return success;
          }

          //
          // Arguments    : void
          // Return Type  : boolean_T
          //
          boolean_T Checkerboard::expandBoardOnce()
          {
            array<double, 3U> b_this;
            array<double, 3U> c_this;
            array<double, 3U> r;
            array<double, 2U> b_index;
            array<double, 2U> b_p2;
            array<double, 2U> idx;
            array<double, 2U> newIndices;
            array<double, 2U> newIndicesLinear;
            array<double, 2U> p1;
            array<double, 2U> p2;
            array<double, 1U> d_this;
            array<double, 1U> e_this;
            array<double, 1U> validIdx;
            array<int, 1U> r1;
            double currCurve_data[5];
            double coordDist;
            double moveDistMultiplier;
            double refCoordValue;
            int i;
            boolean_T success;
            PreviousEnergy = Energy;
            i = 0;
            int exitg1;
            do {
              exitg1 = 0;
              if (i < 4) {
                if (!IsDirectionBad[i]) {
                  float oldEnergy;
                  LastExpandDirection = static_cast<double>(i) + 1.0;
                  oldEnergy = (Energy + static_cast<float>(BoardIdx.size(0) *
                    BoardIdx.size(1))) / static_cast<float>(BoardIdx.size(0) *
                    BoardIdx.size(1));
                  switch (i + 1) {
                   case 1:
                    {
                      int loop_ub;
                      if (IsDistortionHigh) {
                        int numCols;
                        boolean_T exitg2;
                        boolean_T y;
                        fitPolynomialIndices(newIndices);
                        y = true;
                        numCols = 1;
                        exitg2 = false;
                        while ((!exitg2) && (numCols <= newIndices.size(1))) {
                          if (newIndices[numCols - 1] == 0.0) {
                            y = false;
                            exitg2 = true;
                          } else {
                            numCols++;
                          }
                        }

                        if (!y) {
                          int b_loop_ub;
                          b_this.set_size(1, BoardCoords.size(1),
                                          BoardCoords.size(2));
                          loop_ub = BoardCoords.size(2);
                          for (int b_i{0}; b_i < loop_ub; b_i++) {
                            b_loop_ub = BoardCoords.size(1);
                            for (int i1{0}; i1 < b_loop_ub; i1++) {
                              b_this[i1 + b_this.size(1) * b_i] = BoardCoords
                                [(BoardCoords.size(0) * i1 + BoardCoords.size(0)
                                  * BoardCoords.size(1) * b_i) + 1];
                            }
                          }

                          squeeze(b_this, p1);
                          b_this.set_size(1, BoardCoords.size(1),
                                          BoardCoords.size(2));
                          loop_ub = BoardCoords.size(2);
                          for (int b_i{0}; b_i < loop_ub; b_i++) {
                            b_loop_ub = BoardCoords.size(1);
                            for (int i1{0}; i1 < b_loop_ub; i1++) {
                              b_this[i1 + b_this.size(1) * b_i] =
                                BoardCoords[BoardCoords.size(0) * i1 +
                                BoardCoords.size(0) * BoardCoords.size(1) * b_i];
                            }
                          }

                          squeeze(b_this, p2);
                          if ((p2.size(0) == p1.size(0)) && (p2.size(1) ==
                               p1.size(1))) {
                            b_p2.set_size(p2.size(0), p2.size(1));
                            loop_ub = p2.size(0) * p2.size(1);
                            for (int b_i{0}; b_i < loop_ub; b_i++) {
                              b_p2[b_i] = (p2[b_i] + p2[b_i]) - p1[b_i];
                            }

                            findClosestIndices(b_p2, newIndicesLinear);
                          } else {
                            binary_expand_op(newIndicesLinear, this, p2, p1);
                          }

                          numCols = newIndices.size(1);
                          for (b_loop_ub = 0; b_loop_ub < numCols; b_loop_ub++)
                          {
                            if (newIndices[b_loop_ub] == 0.0) {
                              newIndices[b_loop_ub] = newIndicesLinear[b_loop_ub];
                            }
                          }
                        }
                      } else {
                        int b_loop_ub;
                        b_this.set_size(1, BoardCoords.size(1), BoardCoords.size
                                        (2));
                        loop_ub = BoardCoords.size(2);
                        for (int b_i{0}; b_i < loop_ub; b_i++) {
                          b_loop_ub = BoardCoords.size(1);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            b_this[i1 + b_this.size(1) * b_i] = BoardCoords
                              [(BoardCoords.size(0) * i1 + BoardCoords.size(0) *
                                BoardCoords.size(1) * b_i) + 1];
                          }
                        }

                        squeeze(b_this, p1);
                        b_this.set_size(1, BoardCoords.size(1), BoardCoords.size
                                        (2));
                        loop_ub = BoardCoords.size(2);
                        for (int b_i{0}; b_i < loop_ub; b_i++) {
                          b_loop_ub = BoardCoords.size(1);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            b_this[i1 + b_this.size(1) * b_i] =
                              BoardCoords[BoardCoords.size(0) * i1 +
                              BoardCoords.size(0) * BoardCoords.size(1) * b_i];
                          }
                        }

                        squeeze(b_this, p2);
                        if ((p2.size(0) == p1.size(0)) && (p2.size(1) == p1.size
                             (1))) {
                          b_p2.set_size(p2.size(0), p2.size(1));
                          loop_ub = p2.size(0) * p2.size(1);
                          for (int b_i{0}; b_i < loop_ub; b_i++) {
                            b_p2[b_i] = (p2[b_i] + p2[b_i]) - p1[b_i];
                          }

                          findClosestIndices(b_p2, newIndices);
                        } else {
                          binary_expand_op(newIndices, this, p2, p1);
                        }
                      }

                      expandBoardUp(newIndices, p1, r);
                      BoardIdx.set_size(p1.size(0), p1.size(1));
                      loop_ub = p1.size(0) * p1.size(1);
                      for (int b_i{0}; b_i < loop_ub; b_i++) {
                        BoardIdx[b_i] = p1[b_i];
                      }

                      BoardCoords.set_size(r.size(0), r.size(1), r.size(2));
                      loop_ub = r.size(0) * r.size(1) * r.size(2);
                      for (int b_i{0}; b_i < loop_ub; b_i++) {
                        BoardCoords[b_i] = r[b_i];
                      }

                      oldEnergy = computeNewEnergyVertical(oldEnergy);
                    }
                    break;

                   case 2:
                    {
                      int b_i;
                      int loop_ub;
                      int numCols;
                      numCols = BoardCoords.size(0);
                      if (numCols < numCols - 2) {
                        idx.set_size(1, 0);
                      } else {
                        idx.set_size(1, 3);
                        for (b_i = 0; b_i < 3; b_i++) {
                          idx[b_i] = numCols - b_i;
                        }
                      }

                      if (IsDistortionHigh) {
                        double coordsToUse[2];
                        int b_loop_ub;
                        int i1;
                        boolean_T exitg2;
                        boolean_T y;
                        findIndependentVar(idx, coordsToUse);
                        newIndices.set_size(1, BoardCoords.size(1));
                        loop_ub = BoardCoords.size(1);
                        for (b_i = 0; b_i < loop_ub; b_i++) {
                          newIndices[b_i] = 0.0;
                        }

                        newIndicesLinear.set_size(1, 0);
                        b_i = newIndices.size(1);
                        for (int j{0}; j < b_i; j++) {
                          validIdx.set_size(BoardCoords.size(0));
                          loop_ub = BoardCoords.size(0);
                          for (i1 = 0; i1 < loop_ub; i1++) {
                            validIdx[i1] = BoardCoords[(i1 + BoardCoords.size(0)
                              * j) + BoardCoords.size(0) * BoardCoords.size(1) *
                              (static_cast<int>(coordsToUse[0]) - 1)];
                          }

                          eml_find(validIdx, r1);
                          validIdx.set_size(r1.size(0));
                          loop_ub = r1.size(0);
                          for (i1 = 0; i1 < loop_ub; i1++) {
                            validIdx[i1] = r1[i1];
                          }

                          if (validIdx.size(0) >= 2) {
                            double currCoord;
                            double currRad;
                            int currCurve_size[2];
                            findSearchParams(idx, validIdx, static_cast<double>
                                             (j) + 1.0, coordsToUse, &coordDist,
                                             &moveDistMultiplier, &refCoordValue);
                            numCols = 0;
                            i1 = validIdx.size(0);
                            for (b_loop_ub = 0; b_loop_ub < i1; b_loop_ub++) {
                              if (validIdx[b_loop_ub] != 0.0) {
                                numCols++;
                              }
                            }

                            d_this.set_size(validIdx.size(0));
                            loop_ub = validIdx.size(0);
                            for (i1 = 0; i1 < loop_ub; i1++) {
                              d_this[i1] = BoardCoords[((static_cast<int>
                                (validIdx[i1]) + BoardCoords.size(0) * j) +
                                BoardCoords.size(0) * BoardCoords.size(1) * (
                                static_cast<int>(coordsToUse[0]) - 1)) - 1];
                            }

                            e_this.set_size(validIdx.size(0));
                            loop_ub = validIdx.size(0);
                            for (i1 = 0; i1 < loop_ub; i1++) {
                              e_this[i1] = BoardCoords[((static_cast<int>
                                (validIdx[i1]) + BoardCoords.size(0) * j) +
                                BoardCoords.size(0) * BoardCoords.size(1) * (
                                static_cast<int>(coordsToUse[1]) - 1)) - 1];
                            }

                            if (numCols > 5) {
                              i1 = 4;
                            } else {
                              i1 = 2;
                            }

                            polyfit(d_this, e_this, static_cast<double>(i1),
                                    currCurve_data, currCurve_size);
                            currRad = coordDist / 4.0;
                            refCoordValue = BoardCoords[((static_cast<int>
                              (refCoordValue) + BoardCoords.size(0) * j) +
                              BoardCoords.size(0) * BoardCoords.size(1) * (
                              static_cast<int>(coordsToUse[0]) - 1)) - 1];
                            currCoord = currRad + refCoordValue;
                            exitg2 = false;
                            while ((!exitg2) && (std::abs(currCoord -
                                     refCoordValue) < moveDistMultiplier * 1.5 *
                                                 std::abs(coordDist))) {
                              double currPt[2];
                              boolean_T exitg3;
                              y = true;
                              b_loop_ub = 0;
                              exitg3 = false;
                              while ((!exitg3) && (b_loop_ub < 2)) {
                                if (!(coordsToUse[b_loop_ub] == static_cast<
                                      double>(b_loop_ub) + 1.0)) {
                                  y = false;
                                  exitg3 = true;
                                } else {
                                  b_loop_ub++;
                                }
                              }

                              if (y) {
                                double b_y;
                                b_y = currCurve_data[0];
                                i1 = currCurve_size[1];
                                for (b_loop_ub = 0; b_loop_ub <= i1 - 2;
                                     b_loop_ub++) {
                                  b_y = currCoord * b_y +
                                    currCurve_data[b_loop_ub + 1];
                                }

                                currPt[0] = currCoord;
                                currPt[1] = b_y;
                              } else {
                                double b_y;
                                b_y = currCurve_data[0];
                                i1 = currCurve_size[1];
                                for (b_loop_ub = 0; b_loop_ub <= i1 - 2;
                                     b_loop_ub++) {
                                  b_y = currCoord * b_y +
                                    currCurve_data[b_loop_ub + 1];
                                }

                                currPt[0] = b_y;
                                currPt[1] = currCoord;
                              }

                              findClosestOnCurve(currPt, std::abs(currRad),
                                                 currCurve_data, currCurve_size,
                                                 coordsToUse, newIndicesLinear,
                                                 b_index);
                              if (b_index.size(1) != 0) {
                                newIndices[j] = b_index[0];
                                i1 = newIndicesLinear.size(1);
                                loop_ub = b_index.size(1);
                                newIndicesLinear.set_size(newIndicesLinear.size
                                  (0), newIndicesLinear.size(1) + b_index.size(1));
                                for (numCols = 0; numCols < loop_ub; numCols++)
                                {
                                  newIndicesLinear[i1 + numCols] =
                                    b_index[numCols];
                                }

                                exitg2 = true;
                              } else {
                                currCoord += currRad;
                              }
                            }
                          }
                        }

                        y = true;
                        numCols = 1;
                        exitg2 = false;
                        while ((!exitg2) && (numCols <= newIndices.size(1))) {
                          if (newIndices[numCols - 1] == 0.0) {
                            y = false;
                            exitg2 = true;
                          } else {
                            numCols++;
                          }
                        }

                        if (!y) {
                          numCols = static_cast<int>(idx[1]);
                          b_this.set_size(1, BoardCoords.size(1),
                                          BoardCoords.size(2));
                          loop_ub = BoardCoords.size(2);
                          for (b_i = 0; b_i < loop_ub; b_i++) {
                            b_loop_ub = BoardCoords.size(1);
                            for (i1 = 0; i1 < b_loop_ub; i1++) {
                              b_this[i1 + b_this.size(1) * b_i] = BoardCoords
                                [((numCols + BoardCoords.size(0) * i1) +
                                  BoardCoords.size(0) * BoardCoords.size(1) *
                                  b_i) - 1];
                            }
                          }

                          squeeze(b_this, p1);
                          numCols = static_cast<int>(idx[0]);
                          b_this.set_size(1, BoardCoords.size(1),
                                          BoardCoords.size(2));
                          loop_ub = BoardCoords.size(2);
                          for (b_i = 0; b_i < loop_ub; b_i++) {
                            b_loop_ub = BoardCoords.size(1);
                            for (i1 = 0; i1 < b_loop_ub; i1++) {
                              b_this[i1 + b_this.size(1) * b_i] = BoardCoords
                                [((numCols + BoardCoords.size(0) * i1) +
                                  BoardCoords.size(0) * BoardCoords.size(1) *
                                  b_i) - 1];
                            }
                          }

                          squeeze(b_this, p2);
                          if ((p2.size(0) == p1.size(0)) && (p2.size(1) ==
                               p1.size(1))) {
                            b_p2.set_size(p2.size(0), p2.size(1));
                            loop_ub = p2.size(0) * p2.size(1);
                            for (b_i = 0; b_i < loop_ub; b_i++) {
                              b_p2[b_i] = (p2[b_i] + p2[b_i]) - p1[b_i];
                            }

                            findClosestIndices(b_p2, newIndicesLinear);
                          } else {
                            binary_expand_op(newIndicesLinear, this, p2, p1);
                          }

                          numCols = newIndices.size(1);
                          for (b_loop_ub = 0; b_loop_ub < numCols; b_loop_ub++)
                          {
                            if (newIndices[b_loop_ub] == 0.0) {
                              newIndices[b_loop_ub] = newIndicesLinear[b_loop_ub];
                            }
                          }
                        }
                      } else {
                        int b_loop_ub;
                        numCols = static_cast<int>(idx[1]);
                        b_this.set_size(1, BoardCoords.size(1), BoardCoords.size
                                        (2));
                        loop_ub = BoardCoords.size(2);
                        for (b_i = 0; b_i < loop_ub; b_i++) {
                          b_loop_ub = BoardCoords.size(1);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            b_this[i1 + b_this.size(1) * b_i] = BoardCoords
                              [((numCols + BoardCoords.size(0) * i1) +
                                BoardCoords.size(0) * BoardCoords.size(1) * b_i)
                              - 1];
                          }
                        }

                        squeeze(b_this, p1);
                        numCols = static_cast<int>(idx[0]);
                        b_this.set_size(1, BoardCoords.size(1), BoardCoords.size
                                        (2));
                        loop_ub = BoardCoords.size(2);
                        for (b_i = 0; b_i < loop_ub; b_i++) {
                          b_loop_ub = BoardCoords.size(1);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            b_this[i1 + b_this.size(1) * b_i] = BoardCoords
                              [((numCols + BoardCoords.size(0) * i1) +
                                BoardCoords.size(0) * BoardCoords.size(1) * b_i)
                              - 1];
                          }
                        }

                        squeeze(b_this, p2);
                        if ((p2.size(0) == p1.size(0)) && (p2.size(1) == p1.size
                             (1))) {
                          b_p2.set_size(p2.size(0), p2.size(1));
                          loop_ub = p2.size(0) * p2.size(1);
                          for (b_i = 0; b_i < loop_ub; b_i++) {
                            b_p2[b_i] = (p2[b_i] + p2[b_i]) - p1[b_i];
                          }

                          findClosestIndices(b_p2, newIndices);
                        } else {
                          binary_expand_op(newIndices, this, p2, p1);
                        }
                      }

                      expandBoardDown(newIndices, p1, r);
                      BoardIdx.set_size(p1.size(0), p1.size(1));
                      loop_ub = p1.size(0) * p1.size(1);
                      for (b_i = 0; b_i < loop_ub; b_i++) {
                        BoardIdx[b_i] = p1[b_i];
                      }

                      BoardCoords.set_size(r.size(0), r.size(1), r.size(2));
                      loop_ub = r.size(0) * r.size(1) * r.size(2);
                      for (b_i = 0; b_i < loop_ub; b_i++) {
                        BoardCoords[b_i] = r[b_i];
                      }

                      idx.set_size(1, idx.size(1));
                      loop_ub = idx.size(1) - 1;
                      for (b_i = 0; b_i <= loop_ub; b_i++) {
                        idx[b_i] = idx[b_i] + 1.0;
                      }

                      oldEnergy = computeNewEnergyVertical(idx, oldEnergy);
                    }
                    break;

                   case 3:
                    {
                      int loop_ub;
                      if (IsDistortionHigh) {
                        int numCols;
                        boolean_T exitg2;
                        boolean_T y;
                        b_fitPolynomialIndices(newIndices);
                        y = true;
                        numCols = 1;
                        exitg2 = false;
                        while ((!exitg2) && (numCols <= newIndices.size(1))) {
                          if (newIndices[numCols - 1] == 0.0) {
                            y = false;
                            exitg2 = true;
                          } else {
                            numCols++;
                          }
                        }

                        if (!y) {
                          int b_loop_ub;
                          c_this.set_size(BoardCoords.size(0), 1,
                                          BoardCoords.size(2));
                          loop_ub = BoardCoords.size(2);
                          for (int b_i{0}; b_i < loop_ub; b_i++) {
                            b_loop_ub = BoardCoords.size(0);
                            for (int i1{0}; i1 < b_loop_ub; i1++) {
                              c_this[i1 + c_this.size(0) * b_i] = BoardCoords
                                [(i1 + BoardCoords.size(0)) + BoardCoords.size(0)
                                * BoardCoords.size(1) * b_i];
                            }
                          }

                          b_squeeze(c_this, p1);
                          c_this.set_size(BoardCoords.size(0), 1,
                                          BoardCoords.size(2));
                          loop_ub = BoardCoords.size(2);
                          for (int b_i{0}; b_i < loop_ub; b_i++) {
                            b_loop_ub = BoardCoords.size(0);
                            for (int i1{0}; i1 < b_loop_ub; i1++) {
                              c_this[i1 + c_this.size(0) * b_i] = BoardCoords[i1
                                + BoardCoords.size(0) * BoardCoords.size(1) *
                                b_i];
                            }
                          }

                          b_squeeze(c_this, p2);
                          if ((p2.size(0) == p1.size(0)) && (p2.size(1) ==
                               p1.size(1))) {
                            b_p2.set_size(p2.size(0), p2.size(1));
                            loop_ub = p2.size(0) * p2.size(1);
                            for (int b_i{0}; b_i < loop_ub; b_i++) {
                              b_p2[b_i] = (p2[b_i] + p2[b_i]) - p1[b_i];
                            }

                            findClosestIndices(b_p2, newIndicesLinear);
                          } else {
                            binary_expand_op(newIndicesLinear, this, p2, p1);
                          }

                          numCols = newIndices.size(1);
                          for (b_loop_ub = 0; b_loop_ub < numCols; b_loop_ub++)
                          {
                            if (newIndices[b_loop_ub] == 0.0) {
                              newIndices[b_loop_ub] = newIndicesLinear[b_loop_ub];
                            }
                          }
                        }
                      } else {
                        int b_loop_ub;
                        c_this.set_size(BoardCoords.size(0), 1, BoardCoords.size
                                        (2));
                        loop_ub = BoardCoords.size(2);
                        for (int b_i{0}; b_i < loop_ub; b_i++) {
                          b_loop_ub = BoardCoords.size(0);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            c_this[i1 + c_this.size(0) * b_i] = BoardCoords[(i1
                              + BoardCoords.size(0)) + BoardCoords.size(0) *
                              BoardCoords.size(1) * b_i];
                          }
                        }

                        b_squeeze(c_this, p1);
                        c_this.set_size(BoardCoords.size(0), 1, BoardCoords.size
                                        (2));
                        loop_ub = BoardCoords.size(2);
                        for (int b_i{0}; b_i < loop_ub; b_i++) {
                          b_loop_ub = BoardCoords.size(0);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            c_this[i1 + c_this.size(0) * b_i] = BoardCoords[i1 +
                              BoardCoords.size(0) * BoardCoords.size(1) * b_i];
                          }
                        }

                        b_squeeze(c_this, p2);
                        if ((p2.size(0) == p1.size(0)) && (p2.size(1) == p1.size
                             (1))) {
                          b_p2.set_size(p2.size(0), p2.size(1));
                          loop_ub = p2.size(0) * p2.size(1);
                          for (int b_i{0}; b_i < loop_ub; b_i++) {
                            b_p2[b_i] = (p2[b_i] + p2[b_i]) - p1[b_i];
                          }

                          findClosestIndices(b_p2, newIndices);
                        } else {
                          binary_expand_op(newIndices, this, p2, p1);
                        }
                      }

                      expandBoardLeft(newIndices, p1, r);
                      BoardIdx.set_size(p1.size(0), p1.size(1));
                      loop_ub = p1.size(0) * p1.size(1);
                      for (int b_i{0}; b_i < loop_ub; b_i++) {
                        BoardIdx[b_i] = p1[b_i];
                      }

                      BoardCoords.set_size(r.size(0), r.size(1), r.size(2));
                      loop_ub = r.size(0) * r.size(1) * r.size(2);
                      for (int b_i{0}; b_i < loop_ub; b_i++) {
                        BoardCoords[b_i] = r[b_i];
                      }

                      oldEnergy = computeNewEnergyHorizontal(oldEnergy);
                    }
                    break;

                   default:
                    {
                      int loop_ub;
                      int numCols;
                      numCols = BoardCoords.size(1);
                      if (numCols < numCols - 2) {
                        idx.set_size(1, 0);
                      } else {
                        idx.set_size(1, 3);
                        for (int b_i{0}; b_i < 3; b_i++) {
                          idx[b_i] = numCols - b_i;
                        }
                      }

                      if (IsDistortionHigh) {
                        boolean_T exitg2;
                        boolean_T y;
                        fitPolynomialIndices(idx, newIndices);
                        y = true;
                        numCols = 1;
                        exitg2 = false;
                        while ((!exitg2) && (numCols <= newIndices.size(1))) {
                          if (newIndices[numCols - 1] == 0.0) {
                            y = false;
                            exitg2 = true;
                          } else {
                            numCols++;
                          }
                        }

                        if (!y) {
                          int b_loop_ub;
                          numCols = static_cast<int>(idx[1]);
                          c_this.set_size(BoardCoords.size(0), 1,
                                          BoardCoords.size(2));
                          loop_ub = BoardCoords.size(2);
                          for (int b_i{0}; b_i < loop_ub; b_i++) {
                            b_loop_ub = BoardCoords.size(0);
                            for (int i1{0}; i1 < b_loop_ub; i1++) {
                              c_this[i1 + c_this.size(0) * b_i] = BoardCoords
                                [(i1 + BoardCoords.size(0) * (numCols - 1)) +
                                BoardCoords.size(0) * BoardCoords.size(1) * b_i];
                            }
                          }

                          b_squeeze(c_this, p1);
                          numCols = static_cast<int>(idx[0]);
                          c_this.set_size(BoardCoords.size(0), 1,
                                          BoardCoords.size(2));
                          loop_ub = BoardCoords.size(2);
                          for (int b_i{0}; b_i < loop_ub; b_i++) {
                            b_loop_ub = BoardCoords.size(0);
                            for (int i1{0}; i1 < b_loop_ub; i1++) {
                              c_this[i1 + c_this.size(0) * b_i] = BoardCoords
                                [(i1 + BoardCoords.size(0) * (numCols - 1)) +
                                BoardCoords.size(0) * BoardCoords.size(1) * b_i];
                            }
                          }

                          b_squeeze(c_this, p2);
                          if ((p2.size(0) == p1.size(0)) && (p2.size(1) ==
                               p1.size(1))) {
                            b_p2.set_size(p2.size(0), p2.size(1));
                            loop_ub = p2.size(0) * p2.size(1);
                            for (int b_i{0}; b_i < loop_ub; b_i++) {
                              b_p2[b_i] = (p2[b_i] + p2[b_i]) - p1[b_i];
                            }

                            findClosestIndices(b_p2, newIndicesLinear);
                          } else {
                            binary_expand_op(newIndicesLinear, this, p2, p1);
                          }

                          numCols = newIndices.size(1);
                          for (b_loop_ub = 0; b_loop_ub < numCols; b_loop_ub++)
                          {
                            if (newIndices[b_loop_ub] == 0.0) {
                              newIndices[b_loop_ub] = newIndicesLinear[b_loop_ub];
                            }
                          }
                        }
                      } else {
                        int b_loop_ub;
                        numCols = static_cast<int>(idx[1]);
                        c_this.set_size(BoardCoords.size(0), 1, BoardCoords.size
                                        (2));
                        loop_ub = BoardCoords.size(2);
                        for (int b_i{0}; b_i < loop_ub; b_i++) {
                          b_loop_ub = BoardCoords.size(0);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            c_this[i1 + c_this.size(0) * b_i] = BoardCoords[(i1
                              + BoardCoords.size(0) * (numCols - 1)) +
                              BoardCoords.size(0) * BoardCoords.size(1) * b_i];
                          }
                        }

                        b_squeeze(c_this, p1);
                        numCols = static_cast<int>(idx[0]);
                        c_this.set_size(BoardCoords.size(0), 1, BoardCoords.size
                                        (2));
                        loop_ub = BoardCoords.size(2);
                        for (int b_i{0}; b_i < loop_ub; b_i++) {
                          b_loop_ub = BoardCoords.size(0);
                          for (int i1{0}; i1 < b_loop_ub; i1++) {
                            c_this[i1 + c_this.size(0) * b_i] = BoardCoords[(i1
                              + BoardCoords.size(0) * (numCols - 1)) +
                              BoardCoords.size(0) * BoardCoords.size(1) * b_i];
                          }
                        }

                        b_squeeze(c_this, p2);
                        if ((p2.size(0) == p1.size(0)) && (p2.size(1) == p1.size
                             (1))) {
                          b_p2.set_size(p2.size(0), p2.size(1));
                          loop_ub = p2.size(0) * p2.size(1);
                          for (int b_i{0}; b_i < loop_ub; b_i++) {
                            b_p2[b_i] = (p2[b_i] + p2[b_i]) - p1[b_i];
                          }

                          findClosestIndices(b_p2, newIndices);
                        } else {
                          binary_expand_op(newIndices, this, p2, p1);
                        }
                      }

                      expandBoardRight(newIndices, p1, r);
                      BoardIdx.set_size(p1.size(0), p1.size(1));
                      loop_ub = p1.size(0) * p1.size(1);
                      for (int b_i{0}; b_i < loop_ub; b_i++) {
                        BoardIdx[b_i] = p1[b_i];
                      }

                      BoardCoords.set_size(r.size(0), r.size(1), r.size(2));
                      loop_ub = r.size(0) * r.size(1) * r.size(2);
                      for (int b_i{0}; b_i < loop_ub; b_i++) {
                        BoardCoords[b_i] = r[b_i];
                      }

                      idx.set_size(1, idx.size(1));
                      loop_ub = idx.size(1) - 1;
                      for (int b_i{0}; b_i <= loop_ub; b_i++) {
                        idx[b_i] = idx[b_i] + 1.0;
                      }

                      oldEnergy = computeNewEnergyHorizontal(idx, oldEnergy);
                    }
                    break;
                  }

                  Energy = oldEnergy;
                  if (Energy < PreviousEnergy) {
                    success = true;
                    exitg1 = 1;
                  } else {
                    undoLastExpansion();
                    IsDirectionBad[i] = true;
                    i++;
                  }
                } else {
                  i++;
                }
              } else {
                success = false;
                exitg1 = 1;
              }
            } while (exitg1 == 0);

            return success;
          }

          //
          // Arguments    : double seedIdx
          //                const ::coder::array<float, 2U> &points
          //                const float v1[2]
          //                const float v2[2]
          // Return Type  : void
          //
          void Checkerboard::initialize(double seedIdx, const ::coder::array<
            float, 2U> &points, const float v1[2], const float v2[2])
          {
            array<double, 2U> b_r;
            array<float, 2U> center;
            array<float, 2U> d;
            array<float, 2U> l;
            array<float, 2U> pointVectors;
            array<float, 2U> r;
            array<float, 2U> u;
            array<float, 1U> euclideanDists;
            array<boolean_T, 1U> x;
            float b_v1[2];
            int acoef;
            int b_acoef;
            int b_k;
            int csz_idx_1;
            int i;
            int varargin_2;
            int varargin_3;
            boolean_T exitg1;
            boolean_T y;
            BoardIdx.set_size(3, 3);
            for (csz_idx_1 = 0; csz_idx_1 < 9; csz_idx_1++) {
              BoardIdx[csz_idx_1] = 0.0;
            }

            IsDirectionBad[0] = false;
            IsDirectionBad[1] = false;
            IsDirectionBad[2] = false;
            IsDirectionBad[3] = false;
            BoardCoords.set_size(3, 3, 2);
            for (csz_idx_1 = 0; csz_idx_1 < 18; csz_idx_1++) {
              BoardCoords[csz_idx_1] = 0.0;
            }

            Points.set_size(points.size(0), 2);
            acoef = points.size(0) * 2;
            if (static_cast<int>(acoef < 3200)) {
              for (i = 0; i < acoef; i++) {
                Points[i] = points[i];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (i = 0; i < acoef; i++) {
                Points[i] = points[i];
              }
            }

            b_acoef = static_cast<int>(seedIdx);
            center.set_size(1, Points.size(1));
            acoef = Points.size(1);
            if (static_cast<int>(acoef < 3200)) {
              for (i = 0; i < acoef; i++) {
                center[i] = Points[(static_cast<int>(seedIdx) + Points.size(0) *
                                    i) - 1];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (i = 0; i < acoef; i++) {
                center[i] = Points[(b_acoef + Points.size(0) * i) - 1];
              }
            }

            BoardIdx[BoardIdx.size(0) + 1] = seedIdx;
            acoef = BoardCoords.size(2);
            if (static_cast<int>(BoardCoords.size(2) < 3200)) {
              for (i = 0; i < acoef; i++) {
                BoardCoords[(BoardCoords.size(0) + BoardCoords.size(0) *
                             BoardCoords.size(1) * i) + 1] = center[i];
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (i = 0; i < acoef; i++) {
                BoardCoords[(BoardCoords.size(0) + BoardCoords.size(0) *
                             BoardCoords.size(1) * i) + 1] = center[i];
              }
            }

            LastExpandDirection = 1.0;
            PreviousEnergy = rtInfF;
            isValid = false;
            b_acoef = center.size(1);
            acoef = Points.size(1);
            if (b_acoef <= acoef) {
              acoef = b_acoef;
            }

            if (center.size(1) == 1) {
              csz_idx_1 = Points.size(1);
            } else if (Points.size(1) == 1) {
              csz_idx_1 = center.size(1);
            } else if (Points.size(1) == center.size(1)) {
              csz_idx_1 = Points.size(1);
            } else {
              csz_idx_1 = acoef;
            }

            b_acoef = center.size(1);
            acoef = Points.size(1);
            if (b_acoef <= acoef) {
              acoef = b_acoef;
            }

            if (center.size(1) == 1) {
              acoef = Points.size(1);
            } else if (Points.size(1) == 1) {
              acoef = center.size(1);
            } else if (Points.size(1) == center.size(1)) {
              acoef = Points.size(1);
            }

            pointVectors.set_size(Points.size(0), acoef);
            if ((Points.size(0) != 0) && (csz_idx_1 != 0)) {
              int bcoef;
              acoef = (Points.size(1) != 1);
              bcoef = (center.size(1) != 1);
              csz_idx_1--;
              b_acoef = (Points.size(0) != 1);

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32) \
 private(b_k,i,varargin_3,varargin_2)

              for (int k = 0; k <= csz_idx_1; k++) {
                varargin_2 = acoef * k;
                varargin_3 = bcoef * k;
                i = pointVectors.size(0) - 1;
                for (b_k = 0; b_k <= i; b_k++) {
                  pointVectors[b_k + pointVectors.size(0) * k] = Points[b_acoef *
                    b_k + Points.size(0) * varargin_2] - center[varargin_3];
                }
              }
            }

            euclideanDists.set_size(pointVectors.size(0));
            b_acoef = pointVectors.size(0);
            if (static_cast<int>(pointVectors.size(0) < 3200)) {
              for (int k{0}; k < b_acoef; k++) {
                euclideanDists[k] = rt_hypotf_snf(pointVectors[k],
                  pointVectors[k + pointVectors.size(0)]);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (int k = 0; k < b_acoef; k++) {
                euclideanDists[k] = rt_hypotf_snf(pointVectors[k],
                  pointVectors[k + pointVectors.size(0)]);
              }
            }

            BoardIdx[BoardIdx.size(0) * 2 + 1] = findNeighbor(pointVectors,
              euclideanDists, v1);
            b_v1[0] = -v1[0];
            b_v1[1] = -v1[1];
            BoardIdx[1] = findNeighbor(pointVectors, euclideanDists, b_v1);
            BoardIdx[BoardIdx.size(0) + 2] = findNeighbor(pointVectors,
              euclideanDists, v2);
            b_v1[0] = -v2[0];
            b_v1[1] = -v2[1];
            BoardIdx[BoardIdx.size(0)] = findNeighbor(pointVectors,
              euclideanDists, b_v1);
            x.set_size(BoardIdx.size(0) * BoardIdx.size(1));
            acoef = BoardIdx.size(0) * BoardIdx.size(1);
            if (static_cast<int>(acoef < 3200)) {
              for (i = 0; i < acoef; i++) {
                x[i] = (BoardIdx[i] < 0.0);
              }
            } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

              for (i = 0; i < acoef; i++) {
                x[i] = (BoardIdx[i] < 0.0);
              }
            }

            y = false;
            b_acoef = 1;
            exitg1 = false;
            while ((!exitg1) && (b_acoef <= x.size(0))) {
              if (x[b_acoef - 1]) {
                y = true;
                exitg1 = true;
              } else {
                b_acoef++;
              }
            }

            if (y) {
              isValid = false;
            } else {
              b_acoef = static_cast<int>(BoardIdx[BoardIdx.size(0) * 2 + 1]);
              r.set_size(1, Points.size(1));
              acoef = Points.size(1);
              if (static_cast<int>(acoef < 3200)) {
                for (i = 0; i < acoef; i++) {
                  r[i] = Points[(b_acoef + Points.size(0) * i) - 1];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (i = 0; i < acoef; i++) {
                  r[i] = Points[(b_acoef + Points.size(0) * i) - 1];
                }
              }

              acoef = BoardCoords.size(2);
              if (static_cast<int>(BoardCoords.size(2) < 3200)) {
                for (i = 0; i < acoef; i++) {
                  BoardCoords[(BoardCoords.size(0) * 2 + BoardCoords.size(0) *
                               BoardCoords.size(1) * i) + 1] = r[i];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (i = 0; i < acoef; i++) {
                  BoardCoords[(BoardCoords.size(0) * 2 + BoardCoords.size(0) *
                               BoardCoords.size(1) * i) + 1] = r[i];
                }
              }

              b_acoef = static_cast<int>(BoardIdx[1]);
              l.set_size(1, Points.size(1));
              acoef = Points.size(1);
              if (static_cast<int>(acoef < 3200)) {
                for (i = 0; i < acoef; i++) {
                  l[i] = Points[(b_acoef + Points.size(0) * i) - 1];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (i = 0; i < acoef; i++) {
                  l[i] = Points[(b_acoef + Points.size(0) * i) - 1];
                }
              }

              acoef = BoardCoords.size(2);
              if (static_cast<int>(BoardCoords.size(2) < 3200)) {
                for (i = 0; i < acoef; i++) {
                  BoardCoords[BoardCoords.size(0) * BoardCoords.size(1) * i + 1]
                    = l[i];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (i = 0; i < acoef; i++) {
                  BoardCoords[BoardCoords.size(0) * BoardCoords.size(1) * i + 1]
                    = l[i];
                }
              }

              b_acoef = static_cast<int>(BoardIdx[BoardIdx.size(0) + 2]);
              d.set_size(1, Points.size(1));
              acoef = Points.size(1);
              if (static_cast<int>(acoef < 3200)) {
                for (i = 0; i < acoef; i++) {
                  d[i] = Points[(b_acoef + Points.size(0) * i) - 1];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (i = 0; i < acoef; i++) {
                  d[i] = Points[(b_acoef + Points.size(0) * i) - 1];
                }
              }

              acoef = BoardCoords.size(2);
              if (static_cast<int>(BoardCoords.size(2) < 3200)) {
                for (i = 0; i < acoef; i++) {
                  BoardCoords[(BoardCoords.size(0) + BoardCoords.size(0) *
                               BoardCoords.size(1) * i) + 2] = d[i];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (i = 0; i < acoef; i++) {
                  BoardCoords[(BoardCoords.size(0) + BoardCoords.size(0) *
                               BoardCoords.size(1) * i) + 2] = d[i];
                }
              }

              b_acoef = static_cast<int>(BoardIdx[BoardIdx.size(0)]);
              u.set_size(1, Points.size(1));
              acoef = Points.size(1);
              if (static_cast<int>(acoef < 3200)) {
                for (i = 0; i < acoef; i++) {
                  u[i] = Points[(b_acoef + Points.size(0) * i) - 1];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (i = 0; i < acoef; i++) {
                  u[i] = Points[(b_acoef + Points.size(0) * i) - 1];
                }
              }

              acoef = BoardCoords.size(2);
              if (static_cast<int>(BoardCoords.size(2) < 3200)) {
                for (i = 0; i < acoef; i++) {
                  BoardCoords[BoardCoords.size(0) + BoardCoords.size(0) *
                    BoardCoords.size(1) * i] = u[i];
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (i = 0; i < acoef; i++) {
                  BoardCoords[BoardCoords.size(0) + BoardCoords.size(0) *
                    BoardCoords.size(1) * i] = u[i];
                }
              }

              if (u.size(1) == center.size(1)) {
                acoef = u.size(1) - 1;
                u.set_size(1, u.size(1));
                b_acoef = u.size(1) - 1;
                if (static_cast<int>(u.size(1) < 3200)) {
                  for (i = 0; i <= acoef; i++) {
                    u[i] = u[i] - center[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i <= b_acoef; i++) {
                    u[i] = u[i] - center[i];
                  }
                }
              } else {
                minus(u, center);
              }

              if (d.size(1) == center.size(1)) {
                acoef = d.size(1) - 1;
                d.set_size(1, d.size(1));
                b_acoef = d.size(1) - 1;
                if (static_cast<int>(d.size(1) < 3200)) {
                  for (i = 0; i <= acoef; i++) {
                    d[i] = d[i] - center[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i <= b_acoef; i++) {
                    d[i] = d[i] - center[i];
                  }
                }
              } else {
                minus(d, center);
              }

              if (r.size(1) == center.size(1)) {
                acoef = r.size(1) - 1;
                r.set_size(1, r.size(1));
                b_acoef = r.size(1) - 1;
                if (static_cast<int>(r.size(1) < 3200)) {
                  for (i = 0; i <= acoef; i++) {
                    r[i] = r[i] - center[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i <= b_acoef; i++) {
                    r[i] = r[i] - center[i];
                  }
                }
              } else {
                minus(r, center);
              }

              if (l.size(1) == center.size(1)) {
                acoef = l.size(1) - 1;
                l.set_size(1, l.size(1));
                b_acoef = l.size(1) - 1;
                if (static_cast<int>(l.size(1) < 3200)) {
                  for (i = 0; i <= acoef; i++) {
                    l[i] = l[i] - center[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i <= b_acoef; i++) {
                    l[i] = l[i] - center[i];
                  }
                }
              } else {
                minus(l, center);
              }

              if (u.size(1) == l.size(1)) {
                center.set_size(1, u.size(1));
                acoef = u.size(1);
                if (static_cast<int>(u.size(1) < 3200)) {
                  for (i = 0; i < acoef; i++) {
                    center[i] = u[i] + l[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i < acoef; i++) {
                    center[i] = u[i] + l[i];
                  }
                }

                BoardIdx[0] = findNeighbor(pointVectors, euclideanDists, center);
              } else {
                d_binary_expand_op(this, pointVectors, euclideanDists, u, l);
              }

              if (d.size(1) == l.size(1)) {
                center.set_size(1, d.size(1));
                acoef = d.size(1);
                if (static_cast<int>(d.size(1) < 3200)) {
                  for (i = 0; i < acoef; i++) {
                    center[i] = d[i] + l[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i < acoef; i++) {
                    center[i] = d[i] + l[i];
                  }
                }

                BoardIdx[2] = findNeighbor(pointVectors, euclideanDists, center);
              } else {
                c_binary_expand_op(this, pointVectors, euclideanDists, d, l);
              }

              if (d.size(1) == r.size(1)) {
                center.set_size(1, d.size(1));
                acoef = d.size(1);
                if (static_cast<int>(d.size(1) < 3200)) {
                  for (i = 0; i < acoef; i++) {
                    center[i] = d[i] + r[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i < acoef; i++) {
                    center[i] = d[i] + r[i];
                  }
                }

                BoardIdx[BoardIdx.size(0) * 2 + 2] = findNeighbor(pointVectors,
                  euclideanDists, center);
              } else {
                b_binary_expand_op(this, pointVectors, euclideanDists, d, r);
              }

              if (u.size(1) == r.size(1)) {
                center.set_size(1, u.size(1));
                acoef = u.size(1);
                if (static_cast<int>(u.size(1) < 3200)) {
                  for (i = 0; i < acoef; i++) {
                    center[i] = u[i] + r[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i < acoef; i++) {
                    center[i] = u[i] + r[i];
                  }
                }

                BoardIdx[BoardIdx.size(0) * 2] = findNeighbor(pointVectors,
                  euclideanDists, center);
              } else {
                binary_expand_op(this, pointVectors, euclideanDists, u, r);
              }

              x.set_size(BoardIdx.size(0) * BoardIdx.size(1));
              acoef = BoardIdx.size(0) * BoardIdx.size(1);
              if (static_cast<int>(acoef < 3200)) {
                for (i = 0; i < acoef; i++) {
                  x[i] = (BoardIdx[i] > 0.0);
                }
              } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                for (i = 0; i < acoef; i++) {
                  x[i] = (BoardIdx[i] > 0.0);
                }
              }

              isValid = true;
              b_acoef = 1;
              exitg1 = false;
              while ((!exitg1) && (b_acoef <= x.size(0))) {
                if (!x[b_acoef - 1]) {
                  isValid = false;
                  exitg1 = true;
                } else {
                  b_acoef++;
                }
              }

              if (isValid) {
                b_acoef = static_cast<int>(BoardIdx[0]);
                b_r.set_size(1, Points.size(1));
                acoef = Points.size(1);
                if (static_cast<int>(acoef < 3200)) {
                  for (i = 0; i < acoef; i++) {
                    b_r[i] = Points[(b_acoef + Points.size(0) * i) - 1];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i < acoef; i++) {
                    b_r[i] = Points[(b_acoef + Points.size(0) * i) - 1];
                  }
                }

                acoef = BoardCoords.size(2);
                if (static_cast<int>(BoardCoords.size(2) < 3200)) {
                  for (i = 0; i < acoef; i++) {
                    BoardCoords[BoardCoords.size(0) * BoardCoords.size(1) * i] =
                      b_r[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i < acoef; i++) {
                    BoardCoords[BoardCoords.size(0) * BoardCoords.size(1) * i] =
                      b_r[i];
                  }
                }

                b_acoef = static_cast<int>(BoardIdx[2]);
                b_r.set_size(1, Points.size(1));
                acoef = Points.size(1);
                if (static_cast<int>(acoef < 3200)) {
                  for (i = 0; i < acoef; i++) {
                    b_r[i] = Points[(b_acoef + Points.size(0) * i) - 1];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i < acoef; i++) {
                    b_r[i] = Points[(b_acoef + Points.size(0) * i) - 1];
                  }
                }

                acoef = BoardCoords.size(2);
                if (static_cast<int>(BoardCoords.size(2) < 3200)) {
                  for (i = 0; i < acoef; i++) {
                    BoardCoords[BoardCoords.size(0) * BoardCoords.size(1) * i +
                      2] = b_r[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i < acoef; i++) {
                    BoardCoords[BoardCoords.size(0) * BoardCoords.size(1) * i +
                      2] = b_r[i];
                  }
                }

                b_acoef = static_cast<int>(BoardIdx[BoardIdx.size(0) * 2 + 2]);
                b_r.set_size(1, Points.size(1));
                acoef = Points.size(1);
                if (static_cast<int>(acoef < 3200)) {
                  for (i = 0; i < acoef; i++) {
                    b_r[i] = Points[(b_acoef + Points.size(0) * i) - 1];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i < acoef; i++) {
                    b_r[i] = Points[(b_acoef + Points.size(0) * i) - 1];
                  }
                }

                acoef = BoardCoords.size(2);
                if (static_cast<int>(BoardCoords.size(2) < 3200)) {
                  for (i = 0; i < acoef; i++) {
                    BoardCoords[(BoardCoords.size(0) * 2 + BoardCoords.size(0) *
                                 BoardCoords.size(1) * i) + 2] = b_r[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i < acoef; i++) {
                    BoardCoords[(BoardCoords.size(0) * 2 + BoardCoords.size(0) *
                                 BoardCoords.size(1) * i) + 2] = b_r[i];
                  }
                }

                b_acoef = static_cast<int>(BoardIdx[BoardIdx.size(0) * 2]);
                b_r.set_size(1, Points.size(1));
                acoef = Points.size(1);
                if (static_cast<int>(acoef < 3200)) {
                  for (i = 0; i < acoef; i++) {
                    b_r[i] = Points[(b_acoef + Points.size(0) * i) - 1];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i < acoef; i++) {
                    b_r[i] = Points[(b_acoef + Points.size(0) * i) - 1];
                  }
                }

                acoef = BoardCoords.size(2);
                if (static_cast<int>(BoardCoords.size(2) < 3200)) {
                  for (i = 0; i < acoef; i++) {
                    BoardCoords[BoardCoords.size(0) * 2 + BoardCoords.size(0) *
                      BoardCoords.size(1) * i] = b_r[i];
                  }
                } else {

#pragma omp parallel for \
 num_threads(32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

                  for (i = 0; i < acoef; i++) {
                    BoardCoords[BoardCoords.size(0) * 2 + BoardCoords.size(0) *
                      BoardCoords.size(1) * i] = b_r[i];
                  }
                }

                Energy = computeInitialEnergy();
                if (IsDistortionHigh) {
                  acoef = -5;
                } else {
                  acoef = -7;
                }

                isValid = (static_cast<double>(Energy) < acoef);
              }
            }
          }
        }
      }
    }
  }
}

//
// File trailer for Checkerboard.cpp
//
// [EOF]
//
