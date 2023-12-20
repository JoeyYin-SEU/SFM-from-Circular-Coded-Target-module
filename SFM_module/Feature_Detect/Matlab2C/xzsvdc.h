//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xzsvdc.h
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

#ifndef XZSVDC_H
#define XZSVDC_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
namespace internal {
namespace reflapack {
void xzsvdc(const ::coder::array<double, 2U> &A, double *U, double S_data[],
            int *S_size, ::coder::array<double, 2U> &V);

}
} // namespace internal
} // namespace coder

#endif
//
// File trailer for xzsvdc.h
//
// [EOF]
//
