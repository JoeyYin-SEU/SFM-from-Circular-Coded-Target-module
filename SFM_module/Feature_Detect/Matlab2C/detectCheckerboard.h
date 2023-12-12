//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: detectCheckerboard.h
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

#ifndef DETECTCHECKERBOARD_H
#define DETECTCHECKERBOARD_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
namespace coder {
namespace vision {
namespace internal {
namespace calibration {
namespace checkerboard {
class Checkerboard;

}
} // namespace calibration
} // namespace internal
} // namespace vision
} // namespace coder

// Function Declarations
namespace coder {
namespace vision {
namespace internal {
namespace calibration {
namespace checkerboard {
Checkerboard *growCheckerboard(const ::coder::array<float, 2U> &points,
                               const ::coder::array<float, 1U> &scores,
                               const ::coder::array<float, 2U> &Ix2,
                               const ::coder::array<float, 2U> &Iy2,
                               const ::coder::array<float, 2U> &Ixy,
                               double theta, boolean_T highDistortion,
                               boolean_T usePartial, Checkerboard *iobj_0);

Checkerboard *orient(Checkerboard *board, const ::coder::array<float, 2U> &b_I);

void toPoints(const Checkerboard *b_this, boolean_T usePartial,
              ::coder::array<double, 2U> &points, double boardSize[2]);

} // namespace checkerboard
} // namespace calibration
} // namespace internal
} // namespace vision
} // namespace coder

#endif
//
// File trailer for detectCheckerboard.h
//
// [EOF]
//
