//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: subPixelLocation.h
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

#ifndef SUBPIXELLOCATION_H
#define SUBPIXELLOCATION_H

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
void subPixelLocation(const ::coder::array<float, 2U> &metric,
                      ::coder::array<double, 2U> &loc);

}
} // namespace calibration
} // namespace internal
} // namespace vision
} // namespace coder
void subPixelLocationImpl_init();

#endif
//
// File trailer for subPixelLocation.h
//
// [EOF]
//
