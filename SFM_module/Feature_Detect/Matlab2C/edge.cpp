//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: edge.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

// Include Files
#include "edge.h"
#include "get_chessborad_pixel_data.h"
#include "get_chessborad_pixel_rtwutil.h"
#include "imfilter.h"
#include "rt_nonfinite.h"
#include "sum.h"
#include "coder_array.h"
#include "libmwcannythresholding_tbb.h"
#include "libmwedgesobelprewitt_tbb.h"
#include "libmwedgethinning_tbb.h"
#include "libmwgetnumcores.h"
#include "libmwimfilter.h"
#include "libmwippfilter.h"
#include "libmwippreconstruct.h"
#include "libmwtbbhist.h"
#include "omp.h"
#include <cmath>
#include <cstring>

// Function Definitions
//
// Arguments    : const ::coder::array<unsigned char, 2U> &varargin_1
//                ::coder::array<boolean_T, 2U> &varargout_1
// Return Type  : void
//
namespace coder {
void b_edge(const ::coder::array<unsigned char, 2U> &varargin_1,
            ::coder::array<boolean_T, 2U> &varargout_1)
{
  static const double b_kernel[13]{0.0020299751839417137,
                                   0.0102182810143517,
                                   0.058116735860084034,
                                   0.19634433732941292,
                                   0.37823877042972087,
                                   0.35505190018248872,
                                   0.0,
                                   -0.35505190018248872,
                                   -0.37823877042972087,
                                   -0.19634433732941292,
                                   -0.058116735860084034,
                                   -0.0102182810143517,
                                   -0.0020299751839417137};
  static const double c_kernel[13]{
      3.4813359214923059E-5, 0.00054457256285105158, 0.0051667606200595222,
      0.029732654490475543,  0.10377716120747747,    0.21969625200024598,
      0.28209557151935094,   0.21969625200024598,    0.10377716120747747,
      0.029732654490475543,  0.0051667606200595222,  0.00054457256285105158,
      3.4813359214923059E-5};
  static const double d_kernel[13]{0.0020299751839417137,
                                   0.0102182810143517,
                                   0.058116735860084034,
                                   0.19634433732941292,
                                   0.37823877042972087,
                                   0.35505190018248872,
                                   0.0,
                                   -0.35505190018248872,
                                   -0.37823877042972087,
                                   -0.19634433732941292,
                                   -0.058116735860084034,
                                   -0.0102182810143517,
                                   -0.0020299751839417137};
  static const double kernel[13]{
      3.4813359214923059E-5, 0.00054457256285105158, 0.0051667606200595222,
      0.029732654490475543,  0.10377716120747747,    0.21969625200024598,
      0.28209557151935094,   0.21969625200024598,    0.10377716120747747,
      0.029732654490475543,  0.0051667606200595222,  0.00054457256285105158,
      3.4813359214923059E-5};
  static const double nonZeroKernel[12]{
      0.0020299751839417137, 0.0102182810143517,   0.058116735860084034,
      0.19634433732941292,   0.37823877042972087,  0.35505190018248872,
      -0.35505190018248872,  -0.37823877042972087, -0.19634433732941292,
      -0.058116735860084034, -0.0102182810143517,  -0.0020299751839417137};
  static const boolean_T b_conn[13]{true, true, true, true, true, true, false,
                                    true, true, true, true, true, true};
  static const boolean_T c_conn[13]{true, true, true, true, true, true, false,
                                    true, true, true, true, true, true};
  array<float, 2U> a;
  array<float, 2U> b_a;
  array<float, 2U> dx;
  array<boolean_T, 2U> marker;
  array<boolean_T, 2U> mask;
  double counts[64];
  double connDimsT[2];
  double outSizeT[2];
  double padSizeT[2];
  double startT[2];
  double highThreshTemp;
  double numCores;
  float magmax;
  float magmaxPrime;
  int b_m;
  int i;
  int idx;
  int loop_ub;
  int m;
  int n;
  signed char b_idx;
  boolean_T conn[13];
  boolean_T rngFlag;
  boolean_T tooBig;
  a.set_size(varargin_1.size(0), varargin_1.size(1));
  loop_ub = varargin_1.size(0) * varargin_1.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (idx = 0; idx < loop_ub; idx++) {
      a[idx] = static_cast<float>(varargin_1[idx]) / 255.0F;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (idx = 0; idx < loop_ub; idx++) {
      a[idx] = static_cast<float>(varargin_1[idx]) / 255.0F;
    }
  }
  if ((a.size(0) == 0) || (a.size(1) == 0)) {
    varargout_1.set_size(a.size(0), a.size(1));
    loop_ub = a.size(0) * a.size(1);
    if (static_cast<int>(loop_ub < 3200)) {
      for (idx = 0; idx < loop_ub; idx++) {
        varargout_1[idx] = false;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (idx = 0; idx < loop_ub; idx++) {
        varargout_1[idx] = false;
      }
    }
  } else {
    loop_ub = a.size(1);
    m = a.size(0);
    b_m = a.size(0);
    n = a.size(1);
    outSizeT[0] = a.size(0);
    startT[0] = 6.0;
    outSizeT[1] = a.size(1);
    startT[1] = 0.0;
    b_padImage(a, startT, b_a);
    tooBig = true;
    if ((a.size(0) <= 65500) || (a.size(1) <= 65500)) {
      tooBig = false;
    }
    dx.set_size(a.size(0), a.size(1));
    if (!tooBig) {
      padSizeT[0] = b_a.size(0);
      startT[0] = 13.0;
      padSizeT[1] = b_a.size(1);
      startT[1] = 1.0;
      ippfilter_real32(&b_a[0], &dx[0], &outSizeT[0], 2.0, &padSizeT[0],
                       &kernel[0], &startT[0], true);
    } else {
      padSizeT[0] = b_a.size(0);
      padSizeT[1] = b_a.size(1);
      for (i = 0; i < 13; i++) {
        conn[i] = true;
      }
      connDimsT[0] = 13.0;
      connDimsT[1] = 1.0;
      imfilter_real32(&b_a[0], &dx[0], 2.0, &outSizeT[0], 2.0, &padSizeT[0],
                      &kernel[0], 13.0, &conn[0], 2.0, &connDimsT[0],
                      &startT[0], 2.0, true, true);
    }
    outSizeT[0] = dx.size(0);
    startT[0] = 0.0;
    outSizeT[1] = dx.size(1);
    startT[1] = 6.0;
    b_padImage(dx, startT, b_a);
    tooBig = true;
    if ((dx.size(0) <= 65500) || (dx.size(1) <= 65500)) {
      tooBig = false;
    }
    dx.set_size(static_cast<int>(outSizeT[0]), static_cast<int>(outSizeT[1]));
    if (!tooBig) {
      padSizeT[0] = b_a.size(0);
      startT[0] = 1.0;
      padSizeT[1] = b_a.size(1);
      startT[1] = 13.0;
      ippfilter_real32(&b_a[0], &dx[0], &outSizeT[0], 2.0, &padSizeT[0],
                       &b_kernel[0], &startT[0], true);
    } else {
      padSizeT[0] = b_a.size(0);
      connDimsT[0] = 1.0;
      padSizeT[1] = b_a.size(1);
      connDimsT[1] = 13.0;
      imfilter_real32(&b_a[0], &dx[0], 2.0, &outSizeT[0], 2.0, &padSizeT[0],
                      &nonZeroKernel[0], 12.0, &b_conn[0], 2.0, &connDimsT[0],
                      &startT[0], 2.0, true, true);
    }
    outSizeT[0] = a.size(0);
    startT[0] = 0.0;
    outSizeT[1] = a.size(1);
    startT[1] = 6.0;
    b_padImage(a, startT, b_a);
    tooBig = true;
    if ((a.size(0) <= 65500) || (a.size(1) <= 65500)) {
      tooBig = false;
    }
    a.set_size(static_cast<int>(outSizeT[0]), static_cast<int>(outSizeT[1]));
    if (!tooBig) {
      padSizeT[0] = b_a.size(0);
      startT[0] = 1.0;
      padSizeT[1] = b_a.size(1);
      startT[1] = 13.0;
      ippfilter_real32(&b_a[0], &a[0], &outSizeT[0], 2.0, &padSizeT[0],
                       &c_kernel[0], &startT[0], true);
    } else {
      padSizeT[0] = b_a.size(0);
      padSizeT[1] = b_a.size(1);
      for (i = 0; i < 13; i++) {
        conn[i] = true;
      }
      connDimsT[0] = 1.0;
      connDimsT[1] = 13.0;
      imfilter_real32(&b_a[0], &a[0], 2.0, &outSizeT[0], 2.0, &padSizeT[0],
                      &kernel[0], 13.0, &conn[0], 2.0, &connDimsT[0],
                      &startT[0], 2.0, true, true);
    }
    outSizeT[0] = a.size(0);
    startT[0] = 6.0;
    outSizeT[1] = a.size(1);
    startT[1] = 0.0;
    b_padImage(a, startT, b_a);
    tooBig = true;
    if ((a.size(0) <= 65500) || (a.size(1) <= 65500)) {
      tooBig = false;
    }
    a.set_size(static_cast<int>(outSizeT[0]), static_cast<int>(outSizeT[1]));
    if (!tooBig) {
      padSizeT[0] = b_a.size(0);
      startT[0] = 13.0;
      padSizeT[1] = b_a.size(1);
      startT[1] = 1.0;
      ippfilter_real32(&b_a[0], &a[0], &outSizeT[0], 2.0, &padSizeT[0],
                       &d_kernel[0], &startT[0], true);
    } else {
      padSizeT[0] = b_a.size(0);
      connDimsT[0] = 13.0;
      padSizeT[1] = b_a.size(1);
      connDimsT[1] = 1.0;
      imfilter_real32(&b_a[0], &a[0], 2.0, &outSizeT[0], 2.0, &padSizeT[0],
                      &nonZeroKernel[0], 12.0, &c_conn[0], 2.0, &connDimsT[0],
                      &startT[0], 2.0, true, true);
    }
    b_a.set_size(m, loop_ub);
    b_a[0] = rt_hypotf_snf(dx[0], a[0]);
    magmax = b_a[0];
    loop_ub = b_a.size(0) * b_a.size(1) - 2;
#pragma omp parallel num_threads(32 > omp_get_max_threads()                    \
                                     ? omp_get_max_threads()                   \
                                     : 32) private(magmaxPrime)
    {
      magmaxPrime = rtMinusInfF;
#pragma omp for nowait
      for (idx = 0; idx <= loop_ub; idx++) {
        b_a[idx + 1] = rt_hypotf_snf(dx[idx + 1], a[idx + 1]);
        magmaxPrime = std::fmax(b_a[idx + 1], magmaxPrime);
      }
      omp_set_nest_lock(&get_chessborad_pixel_nestLockGlobal);
      {

        magmax = std::fmax(magmaxPrime, magmax);
      }
      omp_unset_nest_lock(&get_chessborad_pixel_nestLockGlobal);
    }
    if (magmax > 0.0F) {
      loop_ub = b_a.size(0) * b_a.size(1);
      if (static_cast<int>(loop_ub < 3200)) {
        for (idx = 0; idx < loop_ub; idx++) {
          b_a[idx] = b_a[idx] / magmax;
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (idx = 0; idx < loop_ub; idx++) {
          b_a[idx] = b_a[idx] / magmax;
        }
      }
    }
    numCores = 1.0;
    getnumcores(&numCores);
    if ((b_a.size(0) * b_a.size(1) > 500000) && (numCores > 1.0)) {
      tooBig = false;
      rngFlag = false;
      tbbhist_real32_scaled(&b_a[0],
                            static_cast<double>(b_a.size(0) * b_a.size(1)),
                            static_cast<double>(b_a.size(0)),
                            static_cast<double>(b_a.size(0) * b_a.size(1)) /
                                static_cast<double>(b_a.size(0)),
                            &counts[0], 64.0, 1.0, 64.0, &rngFlag, &tooBig);
    } else {
      std::memset(&counts[0], 0, 64U * sizeof(double));
      loop_ub = b_a.size(0) * b_a.size(1);
      for (i = 0; i < loop_ub; i++) {
        if (std::isnan(b_a[i])) {
          magmax = 0.0F;
        } else {
          magmax = b_a[i] * 63.0F + 0.5F;
        }
        if (magmax > 63.0F) {
          counts[63]++;
        } else if (std::isinf(b_a[i])) {
          counts[63]++;
        } else {
          counts[static_cast<int>(magmax)]++;
        }
      }
    }
    numCores = 0.0;
    b_idx = 1;
    while ((!(numCores > 0.7 * static_cast<double>(b_a.size(0)) *
                             static_cast<double>(b_a.size(1)))) &&
           (b_idx <= 64)) {
      numCores += counts[b_idx - 1];
      b_idx = static_cast<signed char>(b_idx + 1);
    }
    highThreshTemp = (static_cast<double>(b_idx) - 1.0) / 64.0;
    if ((b_idx > 64) && (!(numCores > 0.7 * static_cast<double>(b_a.size(0)) *
                                          static_cast<double>(b_a.size(1))))) {
    } else {
      numCores = 0.4 * highThreshTemp;
    }
    varargout_1.set_size(b_m, n);
    loop_ub = b_m * n;
    if (static_cast<int>(loop_ub < 3200)) {
      for (idx = 0; idx < loop_ub; idx++) {
        varargout_1[idx] = false;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (idx = 0; idx < loop_ub; idx++) {
        varargout_1[idx] = false;
      }
    }
    if ((b_m != 1) && (n != 1)) {
      startT[0] = b_m;
      startT[1] = n;
      cannythresholding_real32_tbb(&dx[0], &a[0], &b_a[0], &startT[0], numCores,
                                   &varargout_1[0]);
      marker.set_size(b_a.size(0), b_a.size(1));
      loop_ub = b_a.size(0) * b_a.size(1);
      if (static_cast<int>(loop_ub < 3200)) {
        for (idx = 0; idx < loop_ub; idx++) {
          marker[idx] = (b_a[idx] > highThreshTemp);
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (idx = 0; idx < loop_ub; idx++) {
          marker[idx] = (b_a[idx] > highThreshTemp);
        }
      }
      mask.set_size(varargout_1.size(0), varargout_1.size(1));
      loop_ub = varargout_1.size(0) * varargout_1.size(1);
      if (static_cast<int>(loop_ub < 3200)) {
        for (idx = 0; idx < loop_ub; idx++) {
          mask[idx] = varargout_1[idx];
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (idx = 0; idx < loop_ub; idx++) {
          mask[idx] = varargout_1[idx];
        }
      }
      startT[0] = marker.size(0);
      startT[1] = marker.size(1);
      varargout_1.set_size(marker.size(0), marker.size(1));
      loop_ub = marker.size(0) * marker.size(1);
      if (static_cast<int>(loop_ub < 3200)) {
        for (idx = 0; idx < loop_ub; idx++) {
          varargout_1[idx] = marker[idx];
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (idx = 0; idx < loop_ub; idx++) {
          varargout_1[idx] = marker[idx];
        }
      }
      ippreconstruct_uint8((unsigned char *)&varargout_1[0],
                           (unsigned char *)&mask[0], &startT[0], 2.0);
    }
  }
}

//
// Arguments    : const ::coder::array<unsigned char, 2U> &varargin_1
//                ::coder::array<boolean_T, 2U> &varargout_1
// Return Type  : void
//
void c_edge(const ::coder::array<unsigned char, 2U> &varargin_1,
            ::coder::array<boolean_T, 2U> &varargout_1)
{
  static const double kernel[169]{
      5.5729956697509775E-5,  0.00010051139498851065, 0.00020089068977814211,
      0.000368574223248252,   0.00057333707934464928, 0.0007445097188372835,
      0.00081115777138118,    0.0007445097188372835,  0.00057333707934464928,
      0.000368574223248252,   0.00020089068977814211, 0.00010051139498851065,
      5.5729956697509775E-5,  0.00010051139498851065, 0.00023808910723660778,
      0.00052502135634706088, 0.00096021852685293692, 0.0014305914525947186,
      0.0017752325483537438,  0.0018973159045981333,  0.0017752325483537438,
      0.0014305914525947186,  0.00096021852685293692, 0.00052502135634706088,
      0.00023808910723660778, 0.00010051139498851065, 0.00020089068977814211,
      0.00052502135634706088, 0.0011314249555443366,  0.0018973159045981333,
      0.0024906169866851643,  0.0027145805665051842,  0.0027339813609643317,
      0.0027145805665051842,  0.0024906169866851643,  0.0018973159045981333,
      0.0011314249555443366,  0.00052502135634706088, 0.00020089068977814211,
      0.000368574223248252,   0.00096021852685293692, 0.0018973159045981333,
      0.0026624887779519136,  0.0024894667083042682,  0.0014639294238749598,
      0.00084504898031200479, 0.0014639294238749598,  0.0024894667083042682,
      0.0026624887779519136,  0.0018973159045981333,  0.00096021852685293692,
      0.000368574223248252,   0.00057333707934464928, 0.0014305914525947186,
      0.0024906169866851643,  0.0024894667083042682,  3.604838529923914E-5,
      -0.0039654010817205812, -0.0060095999794722091, -0.0039654010817205812,
      3.604838529923914E-5,   0.0024894667083042682,  0.0024906169866851643,
      0.0014305914525947186,  0.00057333707934464928, 0.0007445097188372835,
      0.0017752325483537438,  0.0027145805665051842,  0.0014639294238749598,
      -0.0039654010817205812, -0.011608100855785593,  -0.015357592931151054,
      -0.011608100855785593,  -0.0039654010817205812, 0.0014639294238749598,
      0.0027145805665051842,  0.0017752325483537438,  0.0007445097188372835,
      0.00081115777138118,    0.0018973159045981333,  0.0027339813609643317,
      0.00084504898031200479, -0.0060095999794722091, -0.015357592931151054,
      -0.019899129723045025,  -0.015357592931151054,  -0.0060095999794722091,
      0.00084504898031200479, 0.0027339813609643317,  0.0018973159045981333,
      0.00081115777138118,    0.0007445097188372835,  0.0017752325483537438,
      0.0027145805665051842,  0.0014639294238749598,  -0.0039654010817205812,
      -0.011608100855785593,  -0.015357592931151054,  -0.011608100855785593,
      -0.0039654010817205812, 0.0014639294238749598,  0.0027145805665051842,
      0.0017752325483537438,  0.0007445097188372835,  0.00057333707934464928,
      0.0014305914525947186,  0.0024906169866851643,  0.0024894667083042682,
      3.604838529923914E-5,   -0.0039654010817205812, -0.0060095999794722091,
      -0.0039654010817205812, 3.604838529923914E-5,   0.0024894667083042682,
      0.0024906169866851643,  0.0014305914525947186,  0.00057333707934464928,
      0.000368574223248252,   0.00096021852685293692, 0.0018973159045981333,
      0.0026624887779519136,  0.0024894667083042682,  0.0014639294238749598,
      0.00084504898031200479, 0.0014639294238749598,  0.0024894667083042682,
      0.0026624887779519136,  0.0018973159045981333,  0.00096021852685293692,
      0.000368574223248252,   0.00020089068977814211, 0.00052502135634706088,
      0.0011314249555443366,  0.0018973159045981333,  0.0024906169866851643,
      0.0027145805665051842,  0.0027339813609643317,  0.0027145805665051842,
      0.0024906169866851643,  0.0018973159045981333,  0.0011314249555443366,
      0.00052502135634706088, 0.00020089068977814211, 0.00010051139498851065,
      0.00023808910723660778, 0.00052502135634706088, 0.00096021852685293692,
      0.0014305914525947186,  0.0017752325483537438,  0.0018973159045981333,
      0.0017752325483537438,  0.0014305914525947186,  0.00096021852685293692,
      0.00052502135634706088, 0.00023808910723660778, 0.00010051139498851065,
      5.5729956697509775E-5,  0.00010051139498851065, 0.00020089068977814211,
      0.000368574223248252,   0.00057333707934464928, 0.0007445097188372835,
      0.00081115777138118,    0.0007445097188372835,  0.00057333707934464928,
      0.000368574223248252,   0.00020089068977814211, 0.00010051139498851065,
      5.5729956697509775E-5};
  static const double nonZeroKernel[169]{
      5.5729956697509775E-5,  0.00010051139498851065, 0.00020089068977814211,
      0.000368574223248252,   0.00057333707934464928, 0.0007445097188372835,
      0.00081115777138118,    0.0007445097188372835,  0.00057333707934464928,
      0.000368574223248252,   0.00020089068977814211, 0.00010051139498851065,
      5.5729956697509775E-5,  0.00010051139498851065, 0.00023808910723660778,
      0.00052502135634706088, 0.00096021852685293692, 0.0014305914525947186,
      0.0017752325483537438,  0.0018973159045981333,  0.0017752325483537438,
      0.0014305914525947186,  0.00096021852685293692, 0.00052502135634706088,
      0.00023808910723660778, 0.00010051139498851065, 0.00020089068977814211,
      0.00052502135634706088, 0.0011314249555443366,  0.0018973159045981333,
      0.0024906169866851643,  0.0027145805665051842,  0.0027339813609643317,
      0.0027145805665051842,  0.0024906169866851643,  0.0018973159045981333,
      0.0011314249555443366,  0.00052502135634706088, 0.00020089068977814211,
      0.000368574223248252,   0.00096021852685293692, 0.0018973159045981333,
      0.0026624887779519136,  0.0024894667083042682,  0.0014639294238749598,
      0.00084504898031200479, 0.0014639294238749598,  0.0024894667083042682,
      0.0026624887779519136,  0.0018973159045981333,  0.00096021852685293692,
      0.000368574223248252,   0.00057333707934464928, 0.0014305914525947186,
      0.0024906169866851643,  0.0024894667083042682,  3.604838529923914E-5,
      -0.0039654010817205812, -0.0060095999794722091, -0.0039654010817205812,
      3.604838529923914E-5,   0.0024894667083042682,  0.0024906169866851643,
      0.0014305914525947186,  0.00057333707934464928, 0.0007445097188372835,
      0.0017752325483537438,  0.0027145805665051842,  0.0014639294238749598,
      -0.0039654010817205812, -0.011608100855785593,  -0.015357592931151054,
      -0.011608100855785593,  -0.0039654010817205812, 0.0014639294238749598,
      0.0027145805665051842,  0.0017752325483537438,  0.0007445097188372835,
      0.00081115777138118,    0.0018973159045981333,  0.0027339813609643317,
      0.00084504898031200479, -0.0060095999794722091, -0.015357592931151054,
      -0.019899129723045025,  -0.015357592931151054,  -0.0060095999794722091,
      0.00084504898031200479, 0.0027339813609643317,  0.0018973159045981333,
      0.00081115777138118,    0.0007445097188372835,  0.0017752325483537438,
      0.0027145805665051842,  0.0014639294238749598,  -0.0039654010817205812,
      -0.011608100855785593,  -0.015357592931151054,  -0.011608100855785593,
      -0.0039654010817205812, 0.0014639294238749598,  0.0027145805665051842,
      0.0017752325483537438,  0.0007445097188372835,  0.00057333707934464928,
      0.0014305914525947186,  0.0024906169866851643,  0.0024894667083042682,
      3.604838529923914E-5,   -0.0039654010817205812, -0.0060095999794722091,
      -0.0039654010817205812, 3.604838529923914E-5,   0.0024894667083042682,
      0.0024906169866851643,  0.0014305914525947186,  0.00057333707934464928,
      0.000368574223248252,   0.00096021852685293692, 0.0018973159045981333,
      0.0026624887779519136,  0.0024894667083042682,  0.0014639294238749598,
      0.00084504898031200479, 0.0014639294238749598,  0.0024894667083042682,
      0.0026624887779519136,  0.0018973159045981333,  0.00096021852685293692,
      0.000368574223248252,   0.00020089068977814211, 0.00052502135634706088,
      0.0011314249555443366,  0.0018973159045981333,  0.0024906169866851643,
      0.0027145805665051842,  0.0027339813609643317,  0.0027145805665051842,
      0.0024906169866851643,  0.0018973159045981333,  0.0011314249555443366,
      0.00052502135634706088, 0.00020089068977814211, 0.00010051139498851065,
      0.00023808910723660778, 0.00052502135634706088, 0.00096021852685293692,
      0.0014305914525947186,  0.0017752325483537438,  0.0018973159045981333,
      0.0017752325483537438,  0.0014305914525947186,  0.00096021852685293692,
      0.00052502135634706088, 0.00023808910723660778, 0.00010051139498851065,
      5.5729956697509775E-5,  0.00010051139498851065, 0.00020089068977814211,
      0.000368574223248252,   0.00057333707934464928, 0.0007445097188372835,
      0.00081115777138118,    0.0007445097188372835,  0.00057333707934464928,
      0.000368574223248252,   0.00020089068977814211, 0.00010051139498851065,
      5.5729956697509775E-5};
  array<float, 2U> a;
  array<float, 2U> b_a;
  array<float, 1U> y;
  float b_center;
  float b_down;
  float b_left;
  float b_right;
  float b_up;
  int nx;
  int rr;
  boolean_T zc1;
  boolean_T zc2;
  boolean_T zc3;
  boolean_T zc4;
  a.set_size(varargin_1.size(0), varargin_1.size(1));
  nx = varargin_1.size(0) * varargin_1.size(1);
  if (static_cast<int>(nx < 3200)) {
    for (int k{0}; k < nx; k++) {
      a[k] = static_cast<float>(varargin_1[k]) / 255.0F;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k < nx; k++) {
      a[k] = static_cast<float>(varargin_1[k]) / 255.0F;
    }
  }
  if ((a.size(0) == 0) || (a.size(1) == 0)) {
    varargout_1.set_size(a.size(0), a.size(1));
    nx = a.size(0) * a.size(1);
    if (static_cast<int>(nx < 3200)) {
      for (int k{0}; k < nx; k++) {
        varargout_1[k] = false;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int k = 0; k < nx; k++) {
        varargout_1[k] = false;
      }
    }
  } else {
    double outSizeT[2];
    double startT[2];
    int n;
    int rrmax;
    boolean_T tooBig;
    n = a.size(1);
    rrmax = a.size(0);
    outSizeT[0] = a.size(0);
    startT[0] = 6.0;
    outSizeT[1] = a.size(1);
    startT[1] = 6.0;
    b_padImage(a, startT, b_a);
    tooBig = true;
    if ((a.size(0) <= 65500) || (a.size(1) <= 65500)) {
      tooBig = false;
    }
    a.set_size(static_cast<int>(outSizeT[0]), static_cast<int>(outSizeT[1]));
    if (!tooBig) {
      double padSizeT[2];
      padSizeT[0] = b_a.size(0);
      startT[0] = 13.0;
      padSizeT[1] = b_a.size(1);
      startT[1] = 13.0;
      ippfilter_real32(&b_a[0], &a[0], &outSizeT[0], 2.0, &padSizeT[0],
                       &kernel[0], &startT[0], false);
    } else {
      double padSizeT[2];
      boolean_T conn[169];
      padSizeT[0] = b_a.size(0);
      padSizeT[1] = b_a.size(1);
      for (nx = 0; nx < 169; nx++) {
        conn[nx] = true;
      }
      double connDimsT[2];
      connDimsT[0] = 13.0;
      connDimsT[1] = 13.0;
      imfilter_real32(&b_a[0], &a[0], 2.0, &outSizeT[0], 2.0, &padSizeT[0],
                      &nonZeroKernel[0], 169.0, &conn[0], 2.0, &connDimsT[0],
                      &startT[0], 2.0, true, false);
    }
    nx = a.size(0) * a.size(1);
    y.set_size(a.size(0) * a.size(1));
    if (static_cast<int>(nx < 3200)) {
      for (int k{0}; k < nx; k++) {
        y[k] = std::abs(a[k]);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int k = 0; k < nx; k++) {
        y[k] = std::abs(a[k]);
      }
    }
    startT[0] = 0.75 * sum(y) / static_cast<double>(a.size(0) * a.size(1));
    varargout_1.set_size(rrmax, n);
    nx = rrmax * n;
    if (static_cast<int>(nx < 3200)) {
      for (int k{0}; k < nx; k++) {
        varargout_1[k] = false;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int k = 0; k < nx; k++) {
        varargout_1[k] = false;
      }
    }
    rrmax--;
    nx = n - 3;
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads()                         \
                               : 32) private(zc4, zc3, zc2, zc1, b_center,     \
                                             b_right, b_left, b_down, b_up,    \
                                             rr)

    for (int k = 0; k <= nx; k++) {
      for (rr = 0; rr <= rrmax - 2; rr++) {
        b_up = a[rr + a.size(0) * (k + 1)];
        b_down = a[(rr + a.size(0) * (k + 1)) + 2];
        b_left = a[(rr + a.size(0) * k) + 1];
        b_right = a[(rr + a.size(0) * (k + 2)) + 1];
        b_center = a[(rr + a.size(0) * (k + 1)) + 1];
        if (b_center != 0.0F) {
          if ((b_center < 0.0F) && (b_right > 0.0F) &&
              (std::abs(b_center - b_right) > startT[0])) {
            zc1 = true;
          } else {
            zc1 = false;
          }
          if ((b_left > 0.0F) && (b_center < 0.0F) &&
              (b_left - b_center > startT[0])) {
            zc2 = true;
          } else {
            zc2 = false;
          }
          if ((b_center < 0.0F) && (b_down > 0.0F) &&
              (std::abs(b_center - b_down) > startT[0])) {
            zc3 = true;
          } else {
            zc3 = false;
          }
          if ((b_up > 0.0F) && (b_center < 0.0F) &&
              (b_up - b_center > startT[0])) {
            zc4 = true;
          } else {
            zc4 = false;
          }
        } else {
          if ((b_up < 0.0F) && (b_down > 0.0F) &&
              (std::abs(b_up - b_down) > 2.0 * startT[0])) {
            zc1 = true;
          } else {
            zc1 = false;
          }
          if ((b_up > 0.0F) && (b_down < 0.0F) &&
              (b_up - b_down > 2.0 * startT[0])) {
            zc2 = true;
          } else {
            zc2 = false;
          }
          if ((b_left < 0.0F) && (b_right > 0.0F) &&
              (std::abs(b_left - b_right) > 2.0 * startT[0])) {
            zc3 = true;
          } else {
            zc3 = false;
          }
          if ((b_left > 0.0F) && (b_right < 0.0F) &&
              (b_left - b_right > 2.0 * startT[0])) {
            zc4 = true;
          } else {
            zc4 = false;
          }
        }
        if (zc1 || zc2 || zc3 || zc4) {
          zc4 = true;
        } else {
          zc4 = false;
        }
        varargout_1[(rr + varargout_1.size(0) * (k + 1)) + 1] = zc4;
      }
    }
  }
}

//
// Arguments    : const ::coder::array<unsigned char, 2U> &varargin_1
//                ::coder::array<boolean_T, 2U> &varargout_1
// Return Type  : void
//
void edge(const ::coder::array<unsigned char, 2U> &varargin_1,
          ::coder::array<boolean_T, 2U> &varargout_1)
{
  array<float, 2U> b_b;
  array<float, 2U> pGradient1;
  array<float, 2U> pGradient2;
  array<float, 1U> c_b;
  if ((varargin_1.size(0) == 0) || (varargin_1.size(1) == 0)) {
    int b;
    varargout_1.set_size(varargin_1.size(0), varargin_1.size(1));
    b = varargin_1.size(0) * varargin_1.size(1);
    if (static_cast<int>(b < 3200)) {
      for (int i{0}; i < b; i++) {
        varargout_1[i] = false;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i = 0; i < b; i++) {
        varargout_1[i] = false;
      }
    }
  } else {
    double sz[2];
    int b;
    signed char offset[4];
    sz[0] = varargin_1.size(0);
    sz[1] = varargin_1.size(1);
    pGradient1.set_size(varargin_1.size(0), varargin_1.size(1));
    pGradient2.set_size(varargin_1.size(0), varargin_1.size(1));
    b_b.set_size(varargin_1.size(0), varargin_1.size(1));
    edgesobelprewitt_uint8_tbb(&varargin_1[0], &sz[0], true, 1.0, 1.0,
                               &pGradient1[0], &pGradient2[0], &b_b[0]);
    varargout_1.set_size(varargin_1.size(0), varargin_1.size(1));
    sz[0] = b_b.size(0);
    sz[1] = b_b.size(1);
    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 0;
    offset[3] = 0;
    b = b_b.size(0) * b_b.size(1);
    c_b = b_b.reshape(b);
    edgethinning_real32_tbb(&b_b[0], &pGradient1[0], &pGradient2[0], 1.0, 1.0,
                            &offset[0], 2.2204460492503131E-14,
                            4.0 * sum(c_b) /
                                static_cast<double>(b_b.size(0) * b_b.size(1)),
                            &varargout_1[0], &sz[0]);
  }
}

//
// Arguments    : coder::array<float, 2U> &in1
//                const coder::array<float, 2U> &in2
//                const coder::array<float, 2U> &in3
// Return Type  : void
//
} // namespace coder
void f_binary_expand_op(coder::array<float, 2U> &in1,
                        const coder::array<float, 2U> &in2,
                        const coder::array<float, 2U> &in3)
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
      float f;
      f = in3[i1 * stride_1_0 + in3.size(0) * aux_1_1];
      in1[i1 + in1.size(0) * i] =
          in2[i1 * stride_0_0 + in2.size(0) * aux_0_1] *
              in2[i1 * stride_0_0 + in2.size(0) * aux_0_1] +
          f * f;
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
}

//
// File trailer for edge.cpp
//
// [EOF]
//
