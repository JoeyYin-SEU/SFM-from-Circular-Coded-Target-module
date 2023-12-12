//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: find_peaks.h
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

#ifndef FIND_PEAKS_H
#define FIND_PEAKS_H

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
void find_peaks(const ::coder::array<float, 2U> &metric, double quality,
                ::coder::array<float, 2U> &loc);

}
} // namespace calibration
} // namespace internal
} // namespace vision
} // namespace coder

#endif
//
// File trailer for find_peaks.h
//
// [EOF]
//
