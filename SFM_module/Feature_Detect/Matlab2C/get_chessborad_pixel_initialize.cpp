//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: get_chessborad_pixel_initialize.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "get_chessborad_pixel_initialize.h"
#include "get_chessborad_pixel_data.h"
#include "halideEvalImpl.h"
#include "rt_nonfinite.h"
#include "subPixelLocation.h"
#include "omp.h"

// Function Definitions
//
// Arguments    : void
// Return Type  : void
//
void get_chessborad_pixel_initialize()
{
  omp_init_nest_lock(&get_chessborad_pixel_nestLockGlobal);
  halideInit_not_empty_init();
  subPixelLocationImpl_init();
  isInitialized_get_chessborad_pixel = true;
}

//
// File trailer for get_chessborad_pixel_initialize.cpp
//
// [EOF]
//
