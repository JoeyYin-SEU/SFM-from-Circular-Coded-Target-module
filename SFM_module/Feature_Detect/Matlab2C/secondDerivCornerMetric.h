//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: secondDerivCornerMetric.h
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

#ifndef SECONDDERIVCORNERMETRIC_H
#define SECONDDERIVCORNERMETRIC_H

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
void b_secondDerivCornerMetric(
    const ::coder::array<float, 2U> &b_I, boolean_T highDistortion,
    ::coder::array<float, 2U> &cxy, ::coder::array<float, 2U> &c45,
    ::coder::array<float, 2U> &Ix, ::coder::array<float, 2U> &Iy,
    ::coder::array<float, 2U> &Ixy, ::coder::array<float, 2U> &I_45_45);

void secondDerivCornerMetric(const ::coder::array<float, 2U> &b_I, double sigma,
                             boolean_T highDistortion,
                             ::coder::array<float, 2U> &cxy,
                             ::coder::array<float, 2U> &c45,
                             ::coder::array<float, 2U> &Ix,
                             ::coder::array<float, 2U> &Iy,
                             ::coder::array<float, 2U> &Ixy,
                             ::coder::array<float, 2U> &I_45_45);

} // namespace checkerboard
} // namespace calibration
} // namespace internal
} // namespace vision
} // namespace coder

#endif
//
// File trailer for secondDerivCornerMetric.h
//
// [EOF]
//
