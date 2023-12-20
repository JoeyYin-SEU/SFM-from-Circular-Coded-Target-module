//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: qrsolve.h
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

#ifndef QRSOLVE_H
#define QRSOLVE_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
namespace internal {
void qrsolve(const ::coder::array<double, 2U> &A,
             const ::coder::array<double, 1U> &B, double Y_data[], int *Y_size,
             int *rankA);

}
} // namespace coder

#endif
//
// File trailer for qrsolve.h
//
// [EOF]
//
