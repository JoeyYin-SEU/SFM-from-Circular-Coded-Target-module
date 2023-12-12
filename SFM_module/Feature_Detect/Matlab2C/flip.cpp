//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: flip.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

// Include Files
#include "flip.h"
#include "rt_nonfinite.h"
#include "coder_array.h"

// Function Definitions
//
// Arguments    : ::coder::array<double, 2U> &x
//                double dim
// Return Type  : void
//
namespace coder {
void flip(::coder::array<double, 2U> &x, double dim)
{
  if ((x.size(0) != 0) && (x.size(1) != 0) &&
      (x.size(static_cast<int>(dim) - 1) > 1)) {
    int i;
    int lowerDim;
    int nd2;
    int npages;
    int pagelen;
    int vstride;
    vstride = 1;
    i = static_cast<unsigned char>(static_cast<int>(dim) - 1);
    for (int k{0}; k < i; k++) {
      vstride *= x.size(0);
    }
    pagelen = vstride * x.size(static_cast<int>(dim) - 1);
    npages = 1;
    lowerDim = static_cast<int>(dim) + 1;
    for (int k{lowerDim}; k < 3; k++) {
      npages *= x.size(1);
    }
    i = x.size(static_cast<int>(dim) - 1) - 1;
    nd2 = x.size(static_cast<int>(dim) - 1) >> 1;
    lowerDim = npages - 1;
    for (int j{0}; j <= lowerDim; j++) {
      npages = vstride - 1;
      for (int b_i{0}; b_i <= npages; b_i++) {
        int offset;
        offset = j * pagelen + b_i;
        for (int k{0}; k < nd2; k++) {
          double tmp;
          int i1;
          int tmp_tmp;
          tmp_tmp = offset + k * vstride;
          tmp = x[tmp_tmp];
          i1 = offset + (i - k) * vstride;
          x[tmp_tmp] = x[i1];
          x[i1] = tmp;
        }
      }
    }
  }
}

} // namespace coder

//
// File trailer for flip.cpp
//
// [EOF]
//
