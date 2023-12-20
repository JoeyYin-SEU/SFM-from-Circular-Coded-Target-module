//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: imfilter.h
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

#ifndef IMFILTER_H
#define IMFILTER_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
void b_imfilter(::coder::array<float, 2U> &varargin_1);

void b_padImage(const ::coder::array<float, 2U> &a_tmp, const double pad[2],
                ::coder::array<float, 2U> &a);

void c_imfilter(::coder::array<float, 2U> &varargin_1);

void d_imfilter(::coder::array<float, 2U> &varargin_1);

void imfilter(::coder::array<float, 2U> &varargin_1,
              const double varargin_2_data[], const int varargin_2_size[2]);

void imfilter(::coder::array<float, 2U> &varargin_1,
              const double varargin_2[25]);

void imfilter(::coder::array<float, 2U> &varargin_1);

void imfilter(::coder::array<float, 2U> &varargin_1,
              const ::coder::array<double, 2U> &varargin_2);

void imfilter(::coder::array<float, 2U> &varargin_1,
              const ::coder::array<double, 1U> &varargin_2);

} // namespace coder

#endif
//
// File trailer for imfilter.h
//
// [EOF]
//
