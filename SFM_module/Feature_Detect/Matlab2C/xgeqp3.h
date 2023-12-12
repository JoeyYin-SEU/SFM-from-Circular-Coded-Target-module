//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xgeqp3.h
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

#ifndef XGEQP3_H
#define XGEQP3_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
namespace internal {
namespace lapack {
void xgeqp3(::coder::array<double, 2U> &A, double tau_data[], int *tau_size,
            int jpvt_data[], int jpvt_size[2]);

}
} // namespace internal
} // namespace coder

#endif
//
// File trailer for xgeqp3.h
//
// [EOF]
//
