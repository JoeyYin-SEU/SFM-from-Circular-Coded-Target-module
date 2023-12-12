//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: imedge_2d.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

// Include Files
#include "imedge_2d.h"
#include "edge.h"
#include "get_chessborad_pixel_data.h"
#include "get_chessborad_pixel_initialize.h"
#include "imfilter.h"
#include "rt_nonfinite.h"
#include "sum.h"
#include "coder_array.h"
#include "libmwedgesobelprewitt_tbb.h"
#include "libmwedgethinning_tbb.h"
#include "libmwimfilter.h"
#include "libmwippfilter.h"
#include "omp.h"

// Function Definitions
//
// Arguments    : const coder::array<unsigned char, 2U> &b_I
//                unsigned long long method
//                coder::array<boolean_T, 2U> &bwimage
// Return Type  : void
//
void imedge_2d(const coder::array<unsigned char, 2U> &b_I,
               unsigned long long method, coder::array<boolean_T, 2U> &bwimage)
{
  static const double b_kernel[4]{0.0, -0.5, 0.5, 0.0};
  static const double kernel[4]{0.5, 0.0, 0.0, -0.5};
  static const signed char b_offset[4]{-1, 1, 1, -1};
  static const boolean_T b_conn[4]{false, true, true, false};
  static const boolean_T conn[4]{true, false, false, true};
  coder::array<float, 2U> b_b;
  coder::array<float, 2U> bx;
  coder::array<float, 2U> pGradient1;
  coder::array<float, 2U> pGradient2;
  coder::array<float, 1U> c_b;
  if (!isInitialized_get_chessborad_pixel) {
    get_chessborad_pixel_initialize();
  }
  switch (method) {
  case 0ULL:
    coder::edge(b_I, bwimage);
    break;
  case 1ULL:
    coder::b_edge(b_I, bwimage);
    break;
  case 2ULL: {
    if ((b_I.size(0) == 0) || (b_I.size(1) == 0)) {
      int b;
      bwimage.set_size(b_I.size(0), b_I.size(1));
      b = b_I.size(0) * b_I.size(1);
      if (static_cast<int>(b < 3200)) {
        for (int i{0}; i < b; i++) {
          bwimage[i] = false;
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int i = 0; i < b; i++) {
          bwimage[i] = false;
        }
      }
    } else {
      double sz[2];
      int b;
      signed char offset[4];
      sz[0] = b_I.size(0);
      sz[1] = b_I.size(1);
      pGradient1.set_size(b_I.size(0), b_I.size(1));
      pGradient2.set_size(b_I.size(0), b_I.size(1));
      b_b.set_size(b_I.size(0), b_I.size(1));
      edgesobelprewitt_uint8_tbb(&b_I[0], &sz[0], false, 1.0, 1.0,
                                 &pGradient1[0], &pGradient2[0], &b_b[0]);
      bwimage.set_size(b_I.size(0), b_I.size(1));
      sz[0] = b_b.size(0);
      sz[1] = b_b.size(1);
      offset[0] = 0;
      offset[1] = 0;
      offset[2] = 0;
      offset[3] = 0;
      b = b_b.size(0) * b_b.size(1);
      c_b = b_b.reshape(b);
      edgethinning_real32_tbb(
          &b_b[0], &pGradient1[0], &pGradient2[0], 1.0, 1.0, &offset[0],
          2.2204460492503131E-14,
          4.0 * coder::sum(c_b) /
              static_cast<double>(b_b.size(0) * b_b.size(1)),
          &bwimage[0], &sz[0]);
    }
  } break;
  case 3ULL: {
    int b;
    pGradient2.set_size(b_I.size(0), b_I.size(1));
    b = b_I.size(0) * b_I.size(1);
    if (static_cast<int>(b < 3200)) {
      for (int i{0}; i < b; i++) {
        pGradient2[i] = static_cast<float>(b_I[i]) / 255.0F;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < b; i++) {
        pGradient2[i] = static_cast<float>(b_I[i]) / 255.0F;
      }
    }
    if ((pGradient2.size(0) == 0) || (pGradient2.size(1) == 0)) {
      bwimage.set_size(pGradient2.size(0), pGradient2.size(1));
      b = pGradient2.size(0) * pGradient2.size(1);
      if (static_cast<int>(b < 3200)) {
        for (int i{0}; i < b; i++) {
          bwimage[i] = false;
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int i = 0; i < b; i++) {
          bwimage[i] = false;
        }
      }
    } else {
      double connDimsT[2];
      double nonZeroKernel[2];
      double outSizeT[2];
      double padSizeT[2];
      double sz[2];
      int m;
      int n;
      boolean_T tooBig;
      n = pGradient2.size(1);
      m = pGradient2.size(0);
      outSizeT[0] = pGradient2.size(0);
      sz[0] = 1.0;
      outSizeT[1] = pGradient2.size(1);
      sz[1] = 1.0;
      coder::b_padImage(pGradient2, sz, pGradient1);
      tooBig = true;
      if ((pGradient2.size(0) <= 65500) || (pGradient2.size(1) <= 65500)) {
        tooBig = false;
      }
      bx.set_size(pGradient2.size(0), pGradient2.size(1));
      if (!tooBig) {
        padSizeT[0] = pGradient1.size(0);
        sz[0] = 2.0;
        padSizeT[1] = pGradient1.size(1);
        sz[1] = 2.0;
        ippfilter_real32(&pGradient1[0], &bx[0], &outSizeT[0], 2.0,
                         &padSizeT[0], &kernel[0], &sz[0], false);
      } else {
        padSizeT[0] = pGradient1.size(0);
        nonZeroKernel[0] = 0.5;
        connDimsT[0] = 2.0;
        padSizeT[1] = pGradient1.size(1);
        nonZeroKernel[1] = -0.5;
        connDimsT[1] = 2.0;
        imfilter_real32(&pGradient1[0], &bx[0], 2.0, &outSizeT[0], 2.0,
                        &padSizeT[0], &nonZeroKernel[0], 2.0, &conn[0], 2.0,
                        &connDimsT[0], &sz[0], 2.0, true, false);
      }
      outSizeT[0] = pGradient2.size(0);
      sz[0] = 1.0;
      outSizeT[1] = pGradient2.size(1);
      sz[1] = 1.0;
      coder::b_padImage(pGradient2, sz, pGradient1);
      tooBig = true;
      if ((pGradient2.size(0) <= 65500) || (pGradient2.size(1) <= 65500)) {
        tooBig = false;
      }
      pGradient2.set_size(static_cast<int>(outSizeT[0]),
                          static_cast<int>(outSizeT[1]));
      if (!tooBig) {
        padSizeT[0] = pGradient1.size(0);
        sz[0] = 2.0;
        padSizeT[1] = pGradient1.size(1);
        sz[1] = 2.0;
        ippfilter_real32(&pGradient1[0], &pGradient2[0], &outSizeT[0], 2.0,
                         &padSizeT[0], &b_kernel[0], &sz[0], false);
      } else {
        padSizeT[0] = pGradient1.size(0);
        nonZeroKernel[0] = -0.5;
        connDimsT[0] = 2.0;
        padSizeT[1] = pGradient1.size(1);
        nonZeroKernel[1] = 0.5;
        connDimsT[1] = 2.0;
        imfilter_real32(&pGradient1[0], &pGradient2[0], 2.0, &outSizeT[0], 2.0,
                        &padSizeT[0], &nonZeroKernel[0], 2.0, &b_conn[0], 2.0,
                        &connDimsT[0], &sz[0], 2.0, true, false);
      }
      if ((bx.size(0) == pGradient2.size(0)) &&
          (bx.size(1) == pGradient2.size(1))) {
        b_b.set_size(bx.size(0), bx.size(1));
        b = bx.size(0) * bx.size(1);
        if (static_cast<int>(b < 3200)) {
          for (int i{0}; i < b; i++) {
            b_b[i] = bx[i] * bx[i] + pGradient2[i] * pGradient2[i];
          }
        } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

          for (int i = 0; i < b; i++) {
            b_b[i] = bx[i] * bx[i] + pGradient2[i] * pGradient2[i];
          }
        }
      } else {
        f_binary_expand_op(b_b, bx, pGradient2);
      }
      bwimage.set_size(m, n);
      sz[0] = b_b.size(0);
      sz[1] = b_b.size(1);
      b = b_b.size(0) * b_b.size(1);
      c_b = b_b.reshape(b);
      edgethinning_real32_tbb(
          &b_b[0], &bx[0], &pGradient2[0], 1.0, 1.0, &b_offset[0],
          2.2204460492503131E-14,
          6.0 * coder::sum(c_b) /
              static_cast<double>(b_b.size(0) * b_b.size(1)),
          &bwimage[0], &sz[0]);
    }
  } break;
  case 4ULL:
    coder::c_edge(b_I, bwimage);
    break;
  case 5ULL:
    coder::c_edge(b_I, bwimage);
    break;
  default:
    coder::edge(b_I, bwimage);
    break;
  }
}

//
// File trailer for imedge_2d.cpp
//
// [EOF]
//
