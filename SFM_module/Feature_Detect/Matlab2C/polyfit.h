//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: polyfit.h
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

#ifndef POLYFIT_H
#define POLYFIT_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
void polyfit(const ::coder::array<double, 1U> &x,
             const ::coder::array<double, 1U> &y, double n, double p_data[],
             int p_size[2]);

void polyfit(const ::coder::array<double, 2U> &x,
             const ::coder::array<double, 2U> &y, double n, double p_data[],
             int p_size[2]);

} // namespace coder

#endif
//
// File trailer for polyfit.h
//
// [EOF]
//
