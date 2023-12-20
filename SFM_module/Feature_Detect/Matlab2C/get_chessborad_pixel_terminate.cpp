//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: get_chessborad_pixel_terminate.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "get_chessborad_pixel_terminate.h"
#include "get_chessborad_pixel_data.h"
#include "halideEvalImpl.h"
#include "rt_nonfinite.h"
#include "omp.h"

// Function Declarations
static void customAtExit();

// Function Definitions
//
// Arguments    : void
// Return Type  : void
//
static void customAtExit()
{
  coder::internal::halide::halideCleanup();
}

//
// Arguments    : void
// Return Type  : void
//
void get_chessborad_pixel_terminate()
{
  customAtExit();
  omp_destroy_nest_lock(&get_chessborad_pixel_nestLockGlobal);
  isInitialized_get_chessborad_pixel = false;
}

//
// File trailer for get_chessborad_pixel_terminate.cpp
//
// [EOF]
//
