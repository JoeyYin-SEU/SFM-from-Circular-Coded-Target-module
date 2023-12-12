//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: imdilate.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

// Include Files
#include "imdilate.h"
#include "get_chessborad_pixel_data.h"
#include "rt_nonfinite.h"
#include "HalideRuntime.h"
#include "coder_array.h"
#include "libmwmorphop_ipp.h"
#include "morphop2_halide_imdilate_float32_host.h"
#include "rtGetInf.h"
#include "serializeDeserializeHalideBuffer.hpp"

// Function Definitions
//
// Arguments    : const ::coder::array<float, 2U> &A
//                ::coder::array<float, 2U> &B
// Return Type  : void
//
namespace coder {
void imdilate(const ::coder::array<float, 2U> &A, ::coder::array<float, 2U> &B)
{
  static double minmax[2]{0.0, 0.0};
  minmax[0U] = rtGetMinusInf();
  minmax[1U] = rtGetInf();
  if ((A.size(0) != 0) && (A.size(1) != 0)) {
    halide_buffer_t marshalledInputs[3];
    halide_buffer_t op;
    int ip_size_data[10];
    int size_op[2];
    B.set_size(A.size(0), A.size(1));
    if (!halideInit_not_empty) {
      halideInit_not_empty = true;
      halide_set_num_threads(32);
    }
    ip_size_data[0] = A.size(0);
    ip_size_data[1] = A.size(1);
    marshalledInputs[0] = matlabArrayToHalideBuffer(&A[0], &ip_size_data[0], 2);
    ip_size_data[0] = 2;
    ip_size_data[1] = 1;
    marshalledInputs[1] =
        matlabArrayToHalideBuffer(&minmax[0], &ip_size_data[0], 2);
    ip_size_data[0] = 3;
    ip_size_data[1] = 3;
    marshalledInputs[2] =
        matlabArrayToHalideBuffer(&bv[0], &ip_size_data[0], 2);
    size_op[0] = A.size(0);
    size_op[1] = A.size(1);
    op = matlabArrayToHalideBuffer(&B[0], &size_op[0], 2);
    morphop2_halide_imdilate_float32_host(
        &marshalledInputs[0], &marshalledInputs[1], &marshalledInputs[2], &op);
    halideBufferToMatlabArray(&op, &B[0], B.size(0) * B.size(1));
    deallocateHalideBuffer(&marshalledInputs[0]);
    deallocateHalideBuffer(&marshalledInputs[1]);
    deallocateHalideBuffer(&marshalledInputs[2]);
    deallocateHalideBuffer(&op);
  } else {
    double nsize[2];
    double s[2];
    B.set_size(A.size(0), A.size(1));
    s[0] = A.size(0);
    nsize[0] = 3.0;
    s[1] = A.size(1);
    nsize[1] = 3.0;
    dilate_real32_ipp(&A[0], &s[0], &bv[0], &nsize[0], &B[0]);
  }
}

} // namespace coder

//
// File trailer for imdilate.cpp
//
// [EOF]
//
