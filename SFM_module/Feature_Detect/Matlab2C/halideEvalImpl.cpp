//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: halideEvalImpl.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

// Include Files
#include "halideEvalImpl.h"
#include "get_chessborad_pixel_data.h"
#include "rt_nonfinite.h"
#include "HalideRuntime.h"

// Function Definitions
//
// Arguments    : void
// Return Type  : void
//
namespace coder {
namespace internal {
namespace halide {
void halideCleanup()
{
  halide_shutdown_thread_pool();
}

//
// Arguments    : void
// Return Type  : void
//
} // namespace halide
} // namespace internal
} // namespace coder
void halideInit_not_empty_init()
{
  halideInit_not_empty = false;
}

//
// File trailer for halideEvalImpl.cpp
//
// [EOF]
//
