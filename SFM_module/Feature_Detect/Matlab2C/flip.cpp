//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: flip.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "flip.h"
#include "get_chessborad_pixel_data.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"

// Function Definitions
//
// Arguments    : ::coder::array<double, 2U> &x
//                double dim
// Return Type  : void
//
namespace coder {
void flip(::coder::array<double, 2U> &x, double dim)
{
  double tmp;
  int b_i;
  int i;
  int i1;
  int j;
  int k;
  int lowerDim;
  int nd2;
  int npages;
  int npagesPrime;
  int offset;
  int pagelen;
  int tmp_tmp;
  int vstride;
  if ((x.size(0) != 0) && (x.size(1) != 0) &&
      (x.size(static_cast<int>(dim) - 1) > 1)) {
    vstride = 1;
    i = static_cast<unsigned char>(static_cast<int>(dim) - 1);
    for (k = 0; k < i; k++) {
      vstride *= x.size(0);
    }
    pagelen = vstride * x.size(static_cast<int>(dim) - 1);
    npages = 1;
    lowerDim = static_cast<int>(dim) + 1;
    if (static_cast<int>(2 - static_cast<int>(dim) < 3200)) {
      for (k = lowerDim; k < 3; k++) {
        npages *= x.size(1);
      }
    } else {
#pragma omp parallel num_threads(32 > omp_get_max_threads()                    \
                                     ? omp_get_max_threads()                   \
                                     : 32) private(npagesPrime)
      {
        npagesPrime = 1;
#pragma omp for nowait
        for (k = lowerDim; k < 3; k++) {
          npagesPrime *= x.size(1);
        }
        omp_set_nest_lock(&get_chessborad_pixel_nestLockGlobal);
        {

          npages *= npagesPrime;
        }
        omp_unset_nest_lock(&get_chessborad_pixel_nestLockGlobal);
      }
    }
    i = x.size(static_cast<int>(dim) - 1) - 1;
    nd2 = x.size(static_cast<int>(dim) - 1) >> 1;
    lowerDim = npages - 1;
    for (j = 0; j <= lowerDim; j++) {
      npages = vstride - 1;
      for (b_i = 0; b_i <= npages; b_i++) {
        offset = j * pagelen + b_i;
        for (k = 0; k < nd2; k++) {
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
