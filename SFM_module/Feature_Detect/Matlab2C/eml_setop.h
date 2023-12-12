//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: eml_setop.h
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

#ifndef EML_SETOP_H
#define EML_SETOP_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
void do_vectors(const ::coder::array<double, 2U> &a,
                const ::coder::array<double, 2U> &b,
                ::coder::array<double, 2U> &c, ::coder::array<int, 1U> &ia,
                int *ib_size);

void do_vectors(const ::coder::array<double, 2U> &a,
                const ::coder::array<double, 1U> &b,
                ::coder::array<double, 2U> &c, ::coder::array<int, 1U> &ia,
                int *ib_size);

} // namespace coder

#endif
//
// File trailer for eml_setop.h
//
// [EOF]
//
