//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: get_chessborad_pixel.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

// Include Files
#include "get_chessborad_pixel.h"
#include "Checkerboard.h"
#include "circularBoundary.h"
#include "detectCheckerboard.h"
#include "find_peaks.h"
#include "get_chessborad_pixel_data.h"
#include "get_chessborad_pixel_initialize.h"
#include "imfilter.h"
#include "rt_nonfinite.h"
#include "secondDerivCornerMetric.h"
#include "subPixelLocation.h"
#include "coder_array.h"
#include "omp.h"

// Function Declarations
static void times(coder::array<float, 2U> &in1,
                  const coder::array<float, 2U> &in2);

// Function Definitions
//
// Arguments    : coder::array<float, 2U> &in1
//                const coder::array<float, 2U> &in2
// Return Type  : void
//
static void times(coder::array<float, 2U> &in1,
                  const coder::array<float, 2U> &in2)
{
  coder::array<float, 2U> b_in1;
  int aux_0_1;
  int aux_1_1;
  int c_loop_ub;
  int i;
  int i1;
  int i3;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  if (in2.size(0) == 1) {
    i = in1.size(0);
  } else {
    i = in2.size(0);
  }
  if (in2.size(1) == 1) {
    i1 = in1.size(1);
  } else {
    i1 = in2.size(1);
  }
  b_in1.set_size(i, i1);
  stride_0_0 = (in1.size(0) != 1);
  stride_0_1 = (in1.size(1) != 1);
  stride_1_0 = (in2.size(0) != 1);
  stride_1_1 = (in2.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (in2.size(1) == 1) {
    loop_ub = in1.size(1);
  } else {
    loop_ub = in2.size(1);
  }
  for (i = 0; i < loop_ub; i++) {
    int b_loop_ub;
    i1 = in2.size(0);
    if (i1 == 1) {
      b_loop_ub = in1.size(0);
    } else {
      b_loop_ub = i1;
    }
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_in1[i1 + b_in1.size(0) * i] =
          in1[i1 * stride_0_0 + in1.size(0) * aux_0_1] *
          in2[i1 * stride_1_0 + in2.size(0) * aux_1_1];
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  in1.set_size(b_in1.size(0), b_in1.size(1));
  loop_ub = b_in1.size(1);
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(i3, c_loop_ub)

  for (int i2 = 0; i2 < loop_ub; i2++) {
    c_loop_ub = b_in1.size(0);
    for (i3 = 0; i3 < c_loop_ub; i3++) {
      in1[i3 + in1.size(0) * i2] = b_in1[i3 + b_in1.size(0) * i2];
    }
  }
}

//
// Arguments    : const coder::array<unsigned char, 2U> &b_I
//                double minCornerMetric
//                boolean_T highDistortion
//                boolean_T usePartial
//                coder::array<double, 2U> &imagePoints
//                double boardSize[2]
//                boolean_T *imagesUsed
// Return Type  : void
//
void get_chessborad_pixel(const coder::array<unsigned char, 2U> &b_I,
                          double minCornerMetric, boolean_T highDistortion,
                          boolean_T usePartial,
                          coder::array<double, 2U> &imagePoints,
                          double boardSize[2], boolean_T *imagesUsed)
{
  coder::vision::internal::calibration::checkerboard::Checkerboard b_lobj_0[6];
  coder::vision::internal::calibration::checkerboard::Checkerboard lobj_0[6];
  coder::vision::internal::calibration::checkerboard::Checkerboard *board0;
  coder::vision::internal::calibration::checkerboard::Checkerboard *board45;
  coder::array<double, 1U> s;
  coder::array<float, 2U> I_45_45;
  coder::array<float, 2U> Ix;
  coder::array<float, 2U> Ix2;
  coder::array<float, 2U> Ixy;
  coder::array<float, 2U> Iy;
  coder::array<float, 2U> Iy2;
  coder::array<float, 2U> b_points0;
  coder::array<float, 2U> c45;
  coder::array<float, 2U> c_I;
  coder::array<float, 2U> cxy;
  coder::array<float, 2U> points0;
  coder::array<float, 1U> b_cxy;
  coder::array<int, 1U> r;
  coder::array<boolean_T, 1U> zeroIdx;
  double d;
  float b_varargin_1;
  float varargin_1;
  unsigned int c_varargin_1[2];
  unsigned int varargin_2[2];
  int i1;
  int loop_ub;
  int siz;
  int siz_idx_0;
  boolean_T exitg1;
  boolean_T guard1{false};
  boolean_T guard2{false};
  boolean_T p;
  if (!isInitialized_get_chessborad_pixel) {
    get_chessborad_pixel_initialize();
  }
  c_I.set_size(b_I.size(0), b_I.size(1));
  loop_ub = b_I.size(0) * b_I.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      c_I[i] = static_cast<float>(b_I[i]) / 255.0F;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      c_I[i] = static_cast<float>(b_I[i]) / 255.0F;
    }
  }
  if (highDistortion) {
    d = 1.5;
  } else {
    d = 2.0;
  }
  coder::vision::internal::calibration::checkerboard::secondDerivCornerMetric(
      c_I, d, highDistortion, cxy, c45, Ix, Iy, Ixy, I_45_45);
  Ix2.set_size(Ix.size(0), Ix.size(1));
  loop_ub = Ix.size(0) * Ix.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      varargin_1 = Ix[i];
      Ix2[i] = varargin_1 * varargin_1;
    }
  } else {
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(b_varargin_1)

    for (int i = 0; i < loop_ub; i++) {
      b_varargin_1 = Ix[i];
      Ix2[i] = b_varargin_1 * b_varargin_1;
    }
  }
  coder::c_imfilter(Ix2);
  Iy2.set_size(Iy.size(0), Iy.size(1));
  loop_ub = Iy.size(0) * Iy.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      varargin_1 = Iy[i];
      Iy2[i] = varargin_1 * varargin_1;
    }
  } else {
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(b_varargin_1)

    for (int i = 0; i < loop_ub; i++) {
      b_varargin_1 = Iy[i];
      Iy2[i] = b_varargin_1 * b_varargin_1;
    }
  }
  coder::c_imfilter(Iy2);
  if ((Ix.size(0) == Iy.size(0)) && (Ix.size(1) == Iy.size(1))) {
    loop_ub = Ix.size(0) * Ix.size(1);
    if (static_cast<int>(loop_ub < 3200)) {
      for (int i{0}; i < loop_ub; i++) {
        Ix[i] = Ix[i] * Iy[i];
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < loop_ub; i++) {
        Ix[i] = Ix[i] * Iy[i];
      }
    }
  } else {
    times(Ix, Iy);
  }
  coder::c_imfilter(Ix);
  coder::vision::internal::calibration::checkerboard::find_peaks(
      cxy, minCornerMetric, points0);
  if (highDistortion) {
    b_points0.set_size(points0.size(0), 2);
    loop_ub = points0.size(0) * points0.size(1) - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
      b_points0[i1] = points0[i1];
    }
    coder::vision::internal::calibration::checkerboard::circularBoundary(
        b_points0, c_I, points0);
  }
  siz_idx_0 = cxy.size(0);
  siz = cxy.size(0);
  s.set_size(points0.size(0));
  loop_ub = points0.size(0);
  if (static_cast<int>(points0.size(0) < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      s[i] = static_cast<int>(points0[i + points0.size(0)]) +
             siz_idx_0 * (static_cast<int>(points0[i]) - 1);
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      s[i] = static_cast<int>(points0[i + points0.size(0)]) +
             siz * (static_cast<int>(points0[i]) - 1);
    }
  }
  b_cxy.set_size(s.size(0));
  loop_ub = s.size(0);
  if (static_cast<int>(s.size(0) < 3200)) {
    for (int i{0}; i < loop_ub; i++) {
      b_cxy[i] = cxy[static_cast<int>(s[i]) - 1];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < loop_ub; i++) {
      b_cxy[i] = cxy[static_cast<int>(s[i]) - 1];
    }
  }
  board0 = coder::vision::internal::calibration::checkerboard::growCheckerboard(
      points0, b_cxy, Ix2, Iy2, Ix, 0.0, highDistortion, usePartial,
      &lobj_0[0]);
  imagePoints.set_size(0, 0);
  boardSize[0] = 0.0;
  boardSize[1] = 0.0;
  if (highDistortion) {
    if (board0->isValid) {
      board0 = coder::vision::internal::calibration::checkerboard::orient(
          board0, c_I);
      coder::vision::internal::calibration::checkerboard::toPoints(
          board0, usePartial, imagePoints, boardSize);
      coder::vision::internal::calibration::checkerboard::subPixelLocation(
          Ixy, imagePoints);
    }
  } else {
    coder::vision::internal::calibration::checkerboard::find_peaks(
        c45, minCornerMetric, points0);
    siz_idx_0 = c45.size(0);
    siz = c45.size(0);
    s.set_size(points0.size(0));
    loop_ub = points0.size(0);
    if (static_cast<int>(points0.size(0) < 3200)) {
      for (int i{0}; i < loop_ub; i++) {
        s[i] = static_cast<int>(points0[i + points0.size(0)]) +
               siz_idx_0 * (static_cast<int>(points0[i]) - 1);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < loop_ub; i++) {
        s[i] = static_cast<int>(points0[i + points0.size(0)]) +
               siz * (static_cast<int>(points0[i]) - 1);
      }
    }
    b_cxy.set_size(s.size(0));
    loop_ub = s.size(0);
    if (static_cast<int>(s.size(0) < 3200)) {
      for (int i{0}; i < loop_ub; i++) {
        b_cxy[i] = c45[static_cast<int>(s[i]) - 1];
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < loop_ub; i++) {
        b_cxy[i] = c45[static_cast<int>(s[i]) - 1];
      }
    }
    board45 =
        coder::vision::internal::calibration::checkerboard::growCheckerboard(
            points0, b_cxy, Ix2, Iy2, Ix, 0.78539816339744828, false,
            usePartial, &lobj_0[3]);
    guard1 = false;
    guard2 = false;
    if (board0->isValid) {
      if (board0->Energy <= board45->Energy) {
        guard2 = true;
      } else {
        c_varargin_1[0] = static_cast<unsigned int>(board0->BoardIdx.size(0));
        c_varargin_1[1] = static_cast<unsigned int>(board0->BoardIdx.size(1));
        varargin_2[0] = static_cast<unsigned int>(board45->BoardIdx.size(0));
        varargin_2[1] = static_cast<unsigned int>(board45->BoardIdx.size(1));
        p = true;
        loop_ub = 0;
        exitg1 = false;
        while ((!exitg1) && (loop_ub < 2)) {
          if (static_cast<int>(c_varargin_1[loop_ub]) !=
              static_cast<int>(varargin_2[loop_ub])) {
            p = false;
            exitg1 = true;
          } else {
            loop_ub++;
          }
        }
        if (p) {
          s.set_size(board0->BoardIdx.size(0) * board0->BoardIdx.size(1));
          loop_ub = board0->BoardIdx.size(0) * board0->BoardIdx.size(1);
          for (i1 = 0; i1 < loop_ub; i1++) {
            s[i1] = board0->BoardIdx[i1];
          }
          siz = 0;
          i1 = s.size(0);
          for (loop_ub = 0; loop_ub < i1; loop_ub++) {
            if (s[loop_ub] != 0.0) {
              siz++;
            }
          }
          s.set_size(board45->BoardIdx.size(0) * board45->BoardIdx.size(1));
          loop_ub = board45->BoardIdx.size(0) * board45->BoardIdx.size(1);
          for (i1 = 0; i1 < loop_ub; i1++) {
            s[i1] = board45->BoardIdx[i1];
          }
          siz_idx_0 = 0;
          i1 = s.size(0);
          for (loop_ub = 0; loop_ub < i1; loop_ub++) {
            if (s[loop_ub] != 0.0) {
              siz_idx_0++;
            }
          }
          if (siz > siz_idx_0) {
            guard2 = true;
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }
      }
    } else {
      guard1 = true;
    }
    if (guard2) {
      board0 = coder::vision::internal::calibration::checkerboard::orient(
          board0, c_I);
      coder::vision::internal::calibration::checkerboard::toPoints(
          board0, usePartial, imagePoints, boardSize);
      coder::vision::internal::calibration::checkerboard::subPixelLocation(
          Ixy, imagePoints);
    }
    if (guard1 && board45->isValid) {
      board45 = coder::vision::internal::calibration::checkerboard::orient(
          board45, c_I);
      coder::vision::internal::calibration::checkerboard::toPoints(
          board45, usePartial, imagePoints, boardSize);
      coder::vision::internal::calibration::checkerboard::subPixelLocation(
          I_45_45, imagePoints);
    }
  }
  if ((imagePoints.size(0) == 0) || (imagePoints.size(1) == 0)) {
    coder::vision::internal::calibration::checkerboard::
        b_secondDerivCornerMetric(c_I, highDistortion, cxy, c45, Ix, Iy, Ixy,
                                  I_45_45);
    Ix2.set_size(Ix.size(0), Ix.size(1));
    loop_ub = Ix.size(0) * Ix.size(1);
    if (static_cast<int>(loop_ub < 3200)) {
      for (int i{0}; i < loop_ub; i++) {
        varargin_1 = Ix[i];
        Ix2[i] = varargin_1 * varargin_1;
      }
    } else {
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(b_varargin_1)

      for (int i = 0; i < loop_ub; i++) {
        b_varargin_1 = Ix[i];
        Ix2[i] = b_varargin_1 * b_varargin_1;
      }
    }
    coder::c_imfilter(Ix2);
    Iy2.set_size(Iy.size(0), Iy.size(1));
    loop_ub = Iy.size(0) * Iy.size(1);
    if (static_cast<int>(loop_ub < 3200)) {
      for (int i{0}; i < loop_ub; i++) {
        varargin_1 = Iy[i];
        Iy2[i] = varargin_1 * varargin_1;
      }
    } else {
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(b_varargin_1)

      for (int i = 0; i < loop_ub; i++) {
        b_varargin_1 = Iy[i];
        Iy2[i] = b_varargin_1 * b_varargin_1;
      }
    }
    coder::c_imfilter(Iy2);
    if ((Ix.size(0) == Iy.size(0)) && (Ix.size(1) == Iy.size(1))) {
      loop_ub = Ix.size(0) * Ix.size(1);
      if (static_cast<int>(loop_ub < 3200)) {
        for (int i{0}; i < loop_ub; i++) {
          Ix[i] = Ix[i] * Iy[i];
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int i = 0; i < loop_ub; i++) {
          Ix[i] = Ix[i] * Iy[i];
        }
      }
    } else {
      times(Ix, Iy);
    }
    coder::c_imfilter(Ix);
    coder::vision::internal::calibration::checkerboard::find_peaks(
        cxy, minCornerMetric, points0);
    if (highDistortion) {
      b_points0.set_size(points0.size(0), 2);
      loop_ub = points0.size(0) * points0.size(1) - 1;
      for (i1 = 0; i1 <= loop_ub; i1++) {
        b_points0[i1] = points0[i1];
      }
      coder::vision::internal::calibration::checkerboard::circularBoundary(
          b_points0, c_I, points0);
    }
    siz_idx_0 = cxy.size(0);
    siz = cxy.size(0);
    s.set_size(points0.size(0));
    loop_ub = points0.size(0);
    if (static_cast<int>(points0.size(0) < 3200)) {
      for (int i{0}; i < loop_ub; i++) {
        s[i] = static_cast<int>(points0[i + points0.size(0)]) +
               siz_idx_0 * (static_cast<int>(points0[i]) - 1);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < loop_ub; i++) {
        s[i] = static_cast<int>(points0[i + points0.size(0)]) +
               siz * (static_cast<int>(points0[i]) - 1);
      }
    }
    b_cxy.set_size(s.size(0));
    loop_ub = s.size(0);
    if (static_cast<int>(s.size(0) < 3200)) {
      for (int i{0}; i < loop_ub; i++) {
        b_cxy[i] = cxy[static_cast<int>(s[i]) - 1];
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < loop_ub; i++) {
        b_cxy[i] = cxy[static_cast<int>(s[i]) - 1];
      }
    }
    board0 =
        coder::vision::internal::calibration::checkerboard::growCheckerboard(
            points0, b_cxy, Ix2, Iy2, Ix, 0.0, highDistortion, usePartial,
            &b_lobj_0[0]);
    imagePoints.set_size(0, 0);
    boardSize[0] = 0.0;
    boardSize[1] = 0.0;
    if (highDistortion) {
      if (board0->isValid) {
        board0 = coder::vision::internal::calibration::checkerboard::orient(
            board0, c_I);
        coder::vision::internal::calibration::checkerboard::toPoints(
            board0, usePartial, imagePoints, boardSize);
        coder::vision::internal::calibration::checkerboard::subPixelLocation(
            Ixy, imagePoints);
      }
    } else {
      coder::vision::internal::calibration::checkerboard::find_peaks(
          c45, minCornerMetric, points0);
      siz_idx_0 = c45.size(0);
      siz = c45.size(0);
      s.set_size(points0.size(0));
      loop_ub = points0.size(0);
      if (static_cast<int>(points0.size(0) < 3200)) {
        for (int i{0}; i < loop_ub; i++) {
          s[i] = static_cast<int>(points0[i + points0.size(0)]) +
                 siz_idx_0 * (static_cast<int>(points0[i]) - 1);
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int i = 0; i < loop_ub; i++) {
          s[i] = static_cast<int>(points0[i + points0.size(0)]) +
                 siz * (static_cast<int>(points0[i]) - 1);
        }
      }
      b_cxy.set_size(s.size(0));
      loop_ub = s.size(0);
      if (static_cast<int>(s.size(0) < 3200)) {
        for (int i{0}; i < loop_ub; i++) {
          b_cxy[i] = c45[static_cast<int>(s[i]) - 1];
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int i = 0; i < loop_ub; i++) {
          b_cxy[i] = c45[static_cast<int>(s[i]) - 1];
        }
      }
      board45 =
          coder::vision::internal::calibration::checkerboard::growCheckerboard(
              points0, b_cxy, Ix2, Iy2, Ix, 0.78539816339744828, false,
              usePartial, &b_lobj_0[3]);
      guard1 = false;
      guard2 = false;
      if (board0->isValid) {
        if (board0->Energy <= board45->Energy) {
          guard2 = true;
        } else {
          c_varargin_1[0] = static_cast<unsigned int>(board0->BoardIdx.size(0));
          c_varargin_1[1] = static_cast<unsigned int>(board0->BoardIdx.size(1));
          varargin_2[0] = static_cast<unsigned int>(board45->BoardIdx.size(0));
          varargin_2[1] = static_cast<unsigned int>(board45->BoardIdx.size(1));
          p = true;
          loop_ub = 0;
          exitg1 = false;
          while ((!exitg1) && (loop_ub < 2)) {
            if (static_cast<int>(c_varargin_1[loop_ub]) !=
                static_cast<int>(varargin_2[loop_ub])) {
              p = false;
              exitg1 = true;
            } else {
              loop_ub++;
            }
          }
          if (p) {
            s.set_size(board0->BoardIdx.size(0) * board0->BoardIdx.size(1));
            loop_ub = board0->BoardIdx.size(0) * board0->BoardIdx.size(1);
            for (i1 = 0; i1 < loop_ub; i1++) {
              s[i1] = board0->BoardIdx[i1];
            }
            siz = 0;
            i1 = s.size(0);
            for (loop_ub = 0; loop_ub < i1; loop_ub++) {
              if (s[loop_ub] != 0.0) {
                siz++;
              }
            }
            s.set_size(board45->BoardIdx.size(0) * board45->BoardIdx.size(1));
            loop_ub = board45->BoardIdx.size(0) * board45->BoardIdx.size(1);
            for (i1 = 0; i1 < loop_ub; i1++) {
              s[i1] = board45->BoardIdx[i1];
            }
            siz_idx_0 = 0;
            i1 = s.size(0);
            for (loop_ub = 0; loop_ub < i1; loop_ub++) {
              if (s[loop_ub] != 0.0) {
                siz_idx_0++;
              }
            }
            if (siz > siz_idx_0) {
              guard2 = true;
            } else {
              guard1 = true;
            }
          } else {
            guard1 = true;
          }
        }
      } else {
        guard1 = true;
      }
      if (guard2) {
        board0 = coder::vision::internal::calibration::checkerboard::orient(
            board0, c_I);
        coder::vision::internal::calibration::checkerboard::toPoints(
            board0, usePartial, imagePoints, boardSize);
        coder::vision::internal::calibration::checkerboard::subPixelLocation(
            Ixy, imagePoints);
      }
      if (guard1 && board45->isValid) {
        board45 = coder::vision::internal::calibration::checkerboard::orient(
            board45, c_I);
        coder::vision::internal::calibration::checkerboard::toPoints(
            board45, usePartial, imagePoints, boardSize);
        coder::vision::internal::calibration::checkerboard::subPixelLocation(
            I_45_45, imagePoints);
      }
    }
  }
  if ((imagePoints.size(0) != 0) && (imagePoints.size(1) != 0) && usePartial) {
    zeroIdx.set_size(imagePoints.size(0));
    loop_ub = imagePoints.size(0);
    if (static_cast<int>(imagePoints.size(0) < 3200)) {
      for (int i{0}; i < loop_ub; i++) {
        zeroIdx[i] = (imagePoints[i] == 0.0);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < loop_ub; i++) {
        zeroIdx[i] = (imagePoints[i] == 0.0);
      }
    }
    siz = zeroIdx.size(0) - 1;
    siz_idx_0 = 0;
    for (loop_ub = 0; loop_ub <= siz; loop_ub++) {
      if (zeroIdx[loop_ub]) {
        siz_idx_0++;
      }
    }
    r.set_size(siz_idx_0);
    siz_idx_0 = 0;
    for (loop_ub = 0; loop_ub <= siz; loop_ub++) {
      if (zeroIdx[loop_ub]) {
        r[siz_idx_0] = loop_ub + 1;
        siz_idx_0++;
      }
    }
    loop_ub = imagePoints.size(1);
    siz_idx_0 = r.size(0);
    for (i1 = 0; i1 < loop_ub; i1++) {
      for (siz = 0; siz < siz_idx_0; siz++) {
        imagePoints[(r[siz] + imagePoints.size(0) * i1) - 1] = rtNaN;
      }
    }
  }
  *imagesUsed = ((imagePoints.size(0) != 0) && (imagePoints.size(1) != 0));
}

//
// File trailer for get_chessborad_pixel.cpp
//
// [EOF]
//
