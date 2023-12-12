//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: get_chessborad_pixel.h
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

#ifndef GET_CHESSBORAD_PIXEL_H
#define GET_CHESSBORAD_PIXEL_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
extern void get_chessborad_pixel(const coder::array<unsigned char, 2U> &b_I,
                                 double minCornerMetric,
                                 boolean_T highDistortion, boolean_T usePartial,
                                 coder::array<double, 2U> &imagePoints,
                                 double boardSize[2], boolean_T *imagesUsed);

#endif
//
// File trailer for get_chessborad_pixel.h
//
// [EOF]
//
