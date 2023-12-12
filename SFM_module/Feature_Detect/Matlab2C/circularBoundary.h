//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: circularBoundary.h
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

#ifndef CIRCULARBOUNDARY_H
#define CIRCULARBOUNDARY_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
namespace vision {
namespace internal {
namespace calibration {
namespace checkerboard {
void circularBoundary(const ::coder::array<float, 2U> &imagePoints,
                      const ::coder::array<float, 2U> &b_I,
                      ::coder::array<float, 2U> &points);

}
} // namespace calibration
} // namespace internal
} // namespace vision
} // namespace coder

#endif
//
// File trailer for circularBoundary.h
//
// [EOF]
//
