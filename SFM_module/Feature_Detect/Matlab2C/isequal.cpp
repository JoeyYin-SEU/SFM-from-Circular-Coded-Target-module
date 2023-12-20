//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: isequal.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "isequal.h"
#include "rt_nonfinite.h"

// Function Definitions
//
// Arguments    : const double varargin_1_data[]
//                const int varargin_1_size[2]
//                double varargin_2
// Return Type  : boolean_T
//
namespace coder {
boolean_T isequal(const double varargin_1_data[], const int varargin_1_size[2],
                  double varargin_2)
{
  boolean_T p;
  p = (varargin_1_size[1] == 1);
  if (p && (varargin_1_size[1] != 0) && (!(varargin_1_data[0] == varargin_2))) {
    p = false;
  }
  return p;
}

} // namespace coder

//
// File trailer for isequal.cpp
//
// [EOF]
//
