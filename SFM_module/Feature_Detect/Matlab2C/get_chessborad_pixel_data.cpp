//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: get_chessborad_pixel_data.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "get_chessborad_pixel_data.h"
#include "rt_nonfinite.h"

// Variable Definitions
boolean_T halideInit_not_empty;

omp_nest_lock_t get_chessborad_pixel_nestLockGlobal;

const boolean_T bv[9]{false, true, false, true, true, true, false, true, false};

boolean_T isInitialized_get_chessborad_pixel{false};

//
// File trailer for get_chessborad_pixel_data.cpp
//
// [EOF]
//
