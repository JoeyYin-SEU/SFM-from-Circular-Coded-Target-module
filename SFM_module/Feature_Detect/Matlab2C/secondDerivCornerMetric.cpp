//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: secondDerivCornerMetric.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

// Include Files
#include "secondDerivCornerMetric.h"
#include "bsxfun.h"
#include "colon.h"
#include "combineVectorElements.h"
#include "fspecial.h"
#include "imdilate.h"
#include "imerode.h"
#include "imfilter.h"
#include "minOrMax.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>
#include <cstring>

// Type Definitions
struct cell_wrap_11 {
  double f1[25];
};

// Function Declarations
static void b_binary_expand_op(coder::array<float, 2U> &in1,
                               const coder::array<float, 2U> &in2,
                               const coder::array<float, 2U> &in3);

static void b_binary_expand_op(coder::array<float, 2U> &in1,
                               const coder::array<float, 2U> &in2);

static void binary_expand_op(coder::array<float, 2U> &in1, double in2,
                             double in3, const coder::array<float, 2U> &in4,
                             const coder::array<float, 2U> &in5);

static void binary_expand_op(coder::array<float, 2U> &in1,
                             const coder::array<float, 3U> &in2,
                             const coder::array<float, 2U> &in3);

static void binary_expand_op(coder::array<float, 2U> &in1,
                             const coder::array<float, 3U> &in2);

static void binary_expand_op(coder::array<float, 2U> &in1,
                             const coder::array<float, 2U> &in2,
                             const coder::array<float, 2U> &in3);

static void binary_expand_op(coder::array<float, 2U> &in1,
                             const coder::array<float, 2U> &in2);

static void c_binary_expand_op(coder::array<float, 2U> &in1,
                               const coder::array<float, 2U> &in2);

namespace coder {
namespace vision {
namespace internal {
namespace calibration {
namespace checkerboard {
static void getDerivFilters(cell_wrap_11 filters[12]);

}
} // namespace calibration
} // namespace internal
} // namespace vision
} // namespace coder
static void d_binary_expand_op(coder::array<float, 2U> &in1,
                               const coder::array<float, 2U> &in2,
                               const coder::array<float, 2U> &in3);

// Function Definitions
//
// Arguments    : coder::array<float, 2U> &in1
//                const coder::array<float, 2U> &in2
//                const coder::array<float, 2U> &in3
// Return Type  : void
//
static void b_binary_expand_op(coder::array<float, 2U> &in1,
                               const coder::array<float, 2U> &in2,
                               const coder::array<float, 2U> &in3)
{
  coder::array<float, 2U> b_in1;
  int aux_0_1;
  int aux_1_1;
  int aux_2_1;
  int c_loop_ub;
  int i;
  int i1;
  int i3;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  int stride_2_0;
  int stride_2_1;
  if (in3.size(0) == 1) {
    if (in2.size(0) == 1) {
      i = in1.size(0);
    } else {
      i = in2.size(0);
    }
  } else {
    i = in3.size(0);
  }
  if (in3.size(1) == 1) {
    if (in2.size(1) == 1) {
      i1 = in1.size(1);
    } else {
      i1 = in2.size(1);
    }
  } else {
    i1 = in3.size(1);
  }
  b_in1.set_size(i, i1);
  stride_0_0 = (in1.size(0) != 1);
  stride_0_1 = (in1.size(1) != 1);
  stride_1_0 = (in2.size(0) != 1);
  stride_1_1 = (in2.size(1) != 1);
  stride_2_0 = (in3.size(0) != 1);
  stride_2_1 = (in3.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  aux_2_1 = 0;
  if (in3.size(1) == 1) {
    if (in2.size(1) == 1) {
      loop_ub = in1.size(1);
    } else {
      loop_ub = in2.size(1);
    }
  } else {
    loop_ub = in3.size(1);
  }
  for (i = 0; i < loop_ub; i++) {
    int b_loop_ub;
    i1 = in3.size(0);
    b_loop_ub = in2.size(0);
    if (i1 == 1) {
      if (b_loop_ub == 1) {
        b_loop_ub = in1.size(0);
      }
    } else {
      b_loop_ub = i1;
    }
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_in1[i1 + b_in1.size(0) * i] =
          (in1[i1 * stride_0_0 + in1.size(0) * aux_0_1] +
           in2[i1 * stride_1_0 + in2.size(0) * aux_1_1]) +
          in3[i1 * stride_2_0 + in3.size(0) * aux_2_1];
    }
    aux_2_1 += stride_2_1;
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  in1.set_size(b_in1.size(0), b_in1.size(1));
  loop_ub = b_in1.size(1);
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(i3, c_loop_ub)

  for (int i2 = 0; i2 < loop_ub; i2++) {
    c_loop_ub = b_in1.size(0);
    for (i3 = 0; i3 < c_loop_ub; i3++) {
      in1[i3 + in1.size(0) * i2] = b_in1[i3 + b_in1.size(0) * i2];
    }
  }
}

//
// Arguments    : coder::array<float, 2U> &in1
//                const coder::array<float, 2U> &in2
// Return Type  : void
//
static void b_binary_expand_op(coder::array<float, 2U> &in1,
                               const coder::array<float, 2U> &in2)
{
  coder::array<float, 2U> b_in1;
  int aux_0_1;
  int aux_1_1;
  int c_loop_ub;
  int i;
  int i1;
  int i3;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  if (in2.size(0) == 1) {
    i = in1.size(0);
  } else {
    i = in2.size(0);
  }
  if (in2.size(1) == 1) {
    i1 = in1.size(1);
  } else {
    i1 = in2.size(1);
  }
  b_in1.set_size(i, i1);
  stride_0_0 = (in1.size(0) != 1);
  stride_0_1 = (in1.size(1) != 1);
  stride_1_0 = (in2.size(0) != 1);
  stride_1_1 = (in2.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (in2.size(1) == 1) {
    loop_ub = in1.size(1);
  } else {
    loop_ub = in2.size(1);
  }
  for (i = 0; i < loop_ub; i++) {
    int b_loop_ub;
    i1 = in2.size(0);
    if (i1 == 1) {
      b_loop_ub = in1.size(0);
    } else {
      b_loop_ub = i1;
    }
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_in1[i1 + b_in1.size(0) * i] =
          in1[i1 * stride_0_0 + in1.size(0) * aux_0_1] * 0.707106769F +
          in2[i1 * stride_1_0 + in2.size(0) * aux_1_1] * -0.707106769F;
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  in1.set_size(b_in1.size(0), b_in1.size(1));
  loop_ub = b_in1.size(1);
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(i3, c_loop_ub)

  for (int i2 = 0; i2 < loop_ub; i2++) {
    c_loop_ub = b_in1.size(0);
    for (i3 = 0; i3 < c_loop_ub; i3++) {
      in1[i3 + in1.size(0) * i2] = b_in1[i3 + b_in1.size(0) * i2];
    }
  }
}

//
// Arguments    : coder::array<float, 2U> &in1
//                double in2
//                double in3
//                const coder::array<float, 2U> &in4
//                const coder::array<float, 2U> &in5
// Return Type  : void
//
static void binary_expand_op(coder::array<float, 2U> &in1, double in2,
                             double in3, const coder::array<float, 2U> &in4,
                             const coder::array<float, 2U> &in5)
{
  coder::array<float, 2U> b_in2;
  int aux_0_1;
  int aux_1_1;
  int aux_2_1;
  int c_loop_ub;
  int i;
  int i1;
  int i4;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  int stride_2_0;
  int stride_2_1;
  if (in5.size(0) == 1) {
    i = in4.size(0);
  } else {
    i = in5.size(0);
  }
  if (i == 1) {
    i = in1.size(0);
  } else if (in5.size(0) == 1) {
    i = in4.size(0);
  } else {
    i = in5.size(0);
  }
  if (in5.size(1) == 1) {
    i1 = in4.size(1);
  } else {
    i1 = in5.size(1);
  }
  if (i1 == 1) {
    i1 = in1.size(1);
  } else if (in5.size(1) == 1) {
    i1 = in4.size(1);
  } else {
    i1 = in5.size(1);
  }
  b_in2.set_size(i, i1);
  stride_0_0 = (in1.size(0) != 1);
  stride_0_1 = (in1.size(1) != 1);
  stride_1_0 = (in4.size(0) != 1);
  stride_1_1 = (in4.size(1) != 1);
  stride_2_0 = (in5.size(0) != 1);
  stride_2_1 = (in5.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  aux_2_1 = 0;
  if (in5.size(1) == 1) {
    i = in4.size(1);
  } else {
    i = in5.size(1);
  }
  if (i == 1) {
    loop_ub = in1.size(1);
  } else if (in5.size(1) == 1) {
    loop_ub = in4.size(1);
  } else {
    loop_ub = in5.size(1);
  }
  for (i = 0; i < loop_ub; i++) {
    int b_loop_ub;
    int i3;
    i1 = in5.size(0);
    b_loop_ub = in4.size(0);
    if (i1 == 1) {
      i3 = b_loop_ub;
    } else {
      i3 = i1;
    }
    if (i3 == 1) {
      b_loop_ub = in1.size(0);
    } else if (i1 != 1) {
      b_loop_ub = i1;
    }
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_in2[i1 + b_in2.size(0) * i] =
          static_cast<float>(in2) *
              in1[i1 * stride_0_0 + in1.size(0) * aux_0_1] -
          static_cast<float>(in3) *
              (in4[i1 * stride_1_0 + in4.size(0) * aux_1_1] +
               in5[i1 * stride_2_0 + in5.size(0) * aux_2_1]);
    }
    aux_2_1 += stride_2_1;
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  in1.set_size(b_in2.size(0), b_in2.size(1));
  loop_ub = b_in2.size(1);
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(i4, c_loop_ub)

  for (int i2 = 0; i2 < loop_ub; i2++) {
    c_loop_ub = b_in2.size(0);
    for (i4 = 0; i4 < c_loop_ub; i4++) {
      in1[i4 + in1.size(0) * i2] = b_in2[i4 + b_in2.size(0) * i2];
    }
  }
}

//
// Arguments    : coder::array<float, 2U> &in1
//                const coder::array<float, 3U> &in2
//                const coder::array<float, 2U> &in3
// Return Type  : void
//
static void binary_expand_op(coder::array<float, 2U> &in1,
                             const coder::array<float, 3U> &in2,
                             const coder::array<float, 2U> &in3)
{
  int aux_0_1;
  int aux_1_1;
  int i;
  int i1;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  if (in3.size(0) == 1) {
    i = in2.size(0);
  } else {
    i = in3.size(0);
  }
  if (in3.size(1) == 1) {
    i1 = in2.size(1);
  } else {
    i1 = in3.size(1);
  }
  in1.set_size(i, i1);
  stride_0_0 = (in2.size(0) != 1);
  stride_0_1 = (in2.size(1) != 1);
  stride_1_0 = (in3.size(0) != 1);
  stride_1_1 = (in3.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (in3.size(1) == 1) {
    loop_ub = in2.size(1);
  } else {
    loop_ub = in3.size(1);
  }
  for (i = 0; i < loop_ub; i++) {
    int b_loop_ub;
    i1 = in3.size(0);
    if (i1 == 1) {
      b_loop_ub = in2.size(0);
    } else {
      b_loop_ub = i1;
    }
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      in1[i1 + in1.size(0) * i] =
          (in2[i1 * stride_0_0 + in2.size(0) * aux_0_1] +
           in2[(i1 * stride_0_0 + in2.size(0) * aux_0_1) +
               in2.size(0) * in2.size(1)]) -
          in3[i1 * stride_1_0 + in3.size(0) * aux_1_1];
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
}

//
// Arguments    : coder::array<float, 2U> &in1
//                const coder::array<float, 3U> &in2
// Return Type  : void
//
static void binary_expand_op(coder::array<float, 2U> &in1,
                             const coder::array<float, 3U> &in2)
{
  coder::array<float, 2U> b_in2;
  int aux_0_1;
  int aux_1_1;
  int c_loop_ub;
  int i;
  int i1;
  int i3;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  if (in1.size(0) == 1) {
    i = in2.size(0);
  } else {
    i = in1.size(0);
  }
  if (in1.size(1) == 1) {
    i1 = in2.size(1);
  } else {
    i1 = in1.size(1);
  }
  b_in2.set_size(i, i1);
  stride_0_0 = (in2.size(0) != 1);
  stride_0_1 = (in2.size(1) != 1);
  stride_1_0 = (in1.size(0) != 1);
  stride_1_1 = (in1.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (in1.size(1) == 1) {
    loop_ub = in2.size(1);
  } else {
    loop_ub = in1.size(1);
  }
  for (i = 0; i < loop_ub; i++) {
    int b_loop_ub;
    i1 = in1.size(0);
    if (i1 == 1) {
      b_loop_ub = in2.size(0);
    } else {
      b_loop_ub = i1;
    }
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_in2[i1 + b_in2.size(0) * i] =
          (in2[(i1 * stride_0_0 + in2.size(0) * aux_0_1) +
               in2.size(0) * in2.size(1) * 2] +
           in2[(i1 * stride_0_0 + in2.size(0) * aux_0_1) +
               in2.size(0) * in2.size(1) * 3]) -
          in1[i1 * stride_1_0 + in1.size(0) * aux_1_1];
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  in1.set_size(b_in2.size(0), b_in2.size(1));
  loop_ub = b_in2.size(1);
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(i3, c_loop_ub)

  for (int i2 = 0; i2 < loop_ub; i2++) {
    c_loop_ub = b_in2.size(0);
    for (i3 = 0; i3 < c_loop_ub; i3++) {
      in1[i3 + in1.size(0) * i2] = b_in2[i3 + b_in2.size(0) * i2];
    }
  }
}

//
// Arguments    : coder::array<float, 2U> &in1
//                const coder::array<float, 2U> &in2
//                const coder::array<float, 2U> &in3
// Return Type  : void
//
static void binary_expand_op(coder::array<float, 2U> &in1,
                             const coder::array<float, 2U> &in2,
                             const coder::array<float, 2U> &in3)
{
  int aux_0_1;
  int aux_1_1;
  int i;
  int i1;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  if (in3.size(0) == 1) {
    i = in2.size(0);
  } else {
    i = in3.size(0);
  }
  if (in3.size(1) == 1) {
    i1 = in2.size(1);
  } else {
    i1 = in3.size(1);
  }
  in1.set_size(i, i1);
  stride_0_0 = (in2.size(0) != 1);
  stride_0_1 = (in2.size(1) != 1);
  stride_1_0 = (in3.size(0) != 1);
  stride_1_1 = (in3.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (in3.size(1) == 1) {
    loop_ub = in2.size(1);
  } else {
    loop_ub = in3.size(1);
  }
  for (i = 0; i < loop_ub; i++) {
    int b_loop_ub;
    i1 = in3.size(0);
    if (i1 == 1) {
      b_loop_ub = in2.size(0);
    } else {
      b_loop_ub = i1;
    }
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      in1[i1 + in1.size(0) * i] =
          in2[i1 * stride_0_0 + in2.size(0) * aux_0_1] +
          in3[i1 * stride_1_0 + in3.size(0) * aux_1_1] * 0.707106769F;
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
}

//
// Arguments    : coder::array<float, 2U> &in1
//                const coder::array<float, 2U> &in2
// Return Type  : void
//
static void binary_expand_op(coder::array<float, 2U> &in1,
                             const coder::array<float, 2U> &in2)
{
  coder::array<float, 2U> b_in1;
  int aux_0_1;
  int aux_1_1;
  int c_loop_ub;
  int i;
  int i1;
  int i3;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  if (in2.size(0) == 1) {
    i = in1.size(0);
  } else {
    i = in2.size(0);
  }
  if (in2.size(1) == 1) {
    i1 = in1.size(1);
  } else {
    i1 = in2.size(1);
  }
  b_in1.set_size(i, i1);
  stride_0_0 = (in1.size(0) != 1);
  stride_0_1 = (in1.size(1) != 1);
  stride_1_0 = (in2.size(0) != 1);
  stride_1_1 = (in2.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (in2.size(1) == 1) {
    loop_ub = in1.size(1);
  } else {
    loop_ub = in2.size(1);
  }
  for (i = 0; i < loop_ub; i++) {
    int b_loop_ub;
    i1 = in2.size(0);
    if (i1 == 1) {
      b_loop_ub = in1.size(0);
    } else {
      b_loop_ub = i1;
    }
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_in1[i1 + b_in1.size(0) * i] =
          in1[i1 * stride_0_0 + in1.size(0) * aux_0_1] +
          in2[i1 * stride_1_0 + in2.size(0) * aux_1_1] * -0.707106769F;
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  in1.set_size(b_in1.size(0), b_in1.size(1));
  loop_ub = b_in1.size(1);
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(i3, c_loop_ub)

  for (int i2 = 0; i2 < loop_ub; i2++) {
    c_loop_ub = b_in1.size(0);
    for (i3 = 0; i3 < c_loop_ub; i3++) {
      in1[i3 + in1.size(0) * i2] = b_in1[i3 + b_in1.size(0) * i2];
    }
  }
}

//
// Arguments    : coder::array<float, 2U> &in1
//                const coder::array<float, 2U> &in2
// Return Type  : void
//
static void c_binary_expand_op(coder::array<float, 2U> &in1,
                               const coder::array<float, 2U> &in2)
{
  coder::array<float, 2U> b_in1;
  int aux_0_1;
  int aux_1_1;
  int c_loop_ub;
  int i;
  int i1;
  int i3;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  if (in2.size(0) == 1) {
    i = in1.size(0);
  } else {
    i = in2.size(0);
  }
  if (in2.size(1) == 1) {
    i1 = in1.size(1);
  } else {
    i1 = in2.size(1);
  }
  b_in1.set_size(i, i1);
  stride_0_0 = (in1.size(0) != 1);
  stride_0_1 = (in1.size(1) != 1);
  stride_1_0 = (in2.size(0) != 1);
  stride_1_1 = (in2.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (in2.size(1) == 1) {
    loop_ub = in1.size(1);
  } else {
    loop_ub = in2.size(1);
  }
  for (i = 0; i < loop_ub; i++) {
    int b_loop_ub;
    i1 = in2.size(0);
    if (i1 == 1) {
      b_loop_ub = in1.size(0);
    } else {
      b_loop_ub = i1;
    }
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_in1[i1 + b_in1.size(0) * i] =
          in1[i1 * stride_0_0 + in1.size(0) * aux_0_1] * 0.3F +
          in2[i1 * stride_1_0 + in2.size(0) * aux_1_1] * 0.3F;
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  in1.set_size(b_in1.size(0), b_in1.size(1));
  loop_ub = b_in1.size(1);
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(i3, c_loop_ub)

  for (int i2 = 0; i2 < loop_ub; i2++) {
    c_loop_ub = b_in1.size(0);
    for (i3 = 0; i3 < c_loop_ub; i3++) {
      in1[i3 + in1.size(0) * i2] = b_in1[i3 + b_in1.size(0) * i2];
    }
  }
}

//
// Arguments    : cell_wrap_11 filters[12]
// Return Type  : void
//
namespace coder {
namespace vision {
namespace internal {
namespace calibration {
namespace checkerboard {
static void getDerivFilters(cell_wrap_11 filters[12])
{
  static const double dv[25]{225.0, 206.565051177078,
                             180.0, 153.43494882292202,
                             135.0, 243.43494882292202,
                             225.0, 180.0,
                             135.0, 116.56505117707798,
                             270.0, 270.0,
                             0.0,   90.0,
                             90.0,  296.565051177078,
                             315.0, 0.0,
                             45.0,  63.434948822922024,
                             315.0, 333.434948822922,
                             0.0,   26.565051177077976,
                             45.0};
  static const short angRange[6]{0, 270, 45, 315, 0, 315};
  static const signed char iv[3]{90, 90, 45};
  signed char b_tmp_data[25];
  signed char c_tmp_data[25];
  signed char d_tmp_data[25];
  signed char tmp_data[25];
  std::memset(&filters[0].f1[0], 0, 25U * sizeof(double));
  std::memset(&filters[3].f1[0], 0, 25U * sizeof(double));
  std::memset(&filters[6].f1[0], 0, 25U * sizeof(double));
  std::memset(&filters[9].f1[0], 0, 25U * sizeof(double));
  std::memset(&filters[1].f1[0], 0, 25U * sizeof(double));
  std::memset(&filters[4].f1[0], 0, 25U * sizeof(double));
  std::memset(&filters[7].f1[0], 0, 25U * sizeof(double));
  std::memset(&filters[10].f1[0], 0, 25U * sizeof(double));
  std::memset(&filters[2].f1[0], 0, 25U * sizeof(double));
  std::memset(&filters[5].f1[0], 0, 25U * sizeof(double));
  std::memset(&filters[8].f1[0], 0, 25U * sizeof(double));
  std::memset(&filters[11].f1[0], 0, 25U * sizeof(double));
  for (int rangeIdx{0}; rangeIdx < 3; rangeIdx++) {
    double d;
    double rangeComp;
    double thetaCand_idx_2;
    int partialTrueCount;
    int thetaCand_idx_0;
    int trueCount;
    short i;
    short range_idx_0;
    rangeComp = (180.0 - static_cast<double>(iv[rangeIdx])) / 2.0;
    partialTrueCount = rangeIdx << 1;
    i = angRange[partialTrueCount];
    range_idx_0 = i;
    thetaCand_idx_0 = i;
    i = angRange[partialTrueCount + 1];
    thetaCand_idx_2 = static_cast<double>(range_idx_0) + 2.0 * rangeComp;
    rangeComp = static_cast<double>(i) + -2.0 * rangeComp;
    trueCount = 0;
    partialTrueCount = 0;
    for (int b_i{0}; b_i < 25; b_i++) {
      d = dv[b_i];
      if ((d < thetaCand_idx_0) || (d > i)) {
        trueCount++;
        tmp_data[partialTrueCount] = static_cast<signed char>(b_i + 1);
        partialTrueCount++;
      }
    }
    for (partialTrueCount = 0; partialTrueCount < trueCount;
         partialTrueCount++) {
      filters[rangeIdx].f1[tmp_data[partialTrueCount] - 1] = 1.0;
    }
    trueCount = 0;
    partialTrueCount = 0;
    for (int b_i{0}; b_i < 25; b_i++) {
      d = dv[b_i];
      if ((d > thetaCand_idx_2) && (d < rangeComp)) {
        trueCount++;
        b_tmp_data[partialTrueCount] = static_cast<signed char>(b_i + 1);
        partialTrueCount++;
      }
    }
    for (partialTrueCount = 0; partialTrueCount < trueCount;
         partialTrueCount++) {
      filters[rangeIdx].f1[b_tmp_data[partialTrueCount] - 1] = -1.0;
    }
    filters[rangeIdx].f1[12] = 0.0;
    trueCount = 0;
    partialTrueCount = 0;
    for (int b_i{0}; b_i < 25; b_i++) {
      d = dv[b_i];
      if ((d > rangeComp) && (d < i)) {
        trueCount++;
        c_tmp_data[partialTrueCount] = static_cast<signed char>(b_i + 1);
        partialTrueCount++;
      }
    }
    for (partialTrueCount = 0; partialTrueCount < trueCount;
         partialTrueCount++) {
      filters[rangeIdx + 3].f1[c_tmp_data[partialTrueCount] - 1] = 1.0;
    }
    trueCount = 0;
    partialTrueCount = 0;
    for (int b_i{0}; b_i < 25; b_i++) {
      d = dv[b_i];
      if ((d > thetaCand_idx_0) && (d < thetaCand_idx_2)) {
        trueCount++;
        d_tmp_data[partialTrueCount] = static_cast<signed char>(b_i + 1);
        partialTrueCount++;
      }
    }
    for (partialTrueCount = 0; partialTrueCount < trueCount;
         partialTrueCount++) {
      filters[rangeIdx + 3].f1[d_tmp_data[partialTrueCount] - 1] = -1.0;
    }
    filters[rangeIdx + 3].f1[12] = 0.0;
    for (partialTrueCount = 0; partialTrueCount < 25; partialTrueCount++) {
      filters[rangeIdx + 6].f1[partialTrueCount] =
          -filters[rangeIdx].f1[partialTrueCount];
      filters[rangeIdx + 9].f1[partialTrueCount] =
          -filters[rangeIdx + 3].f1[partialTrueCount];
    }
  }
}

//
// Arguments    : coder::array<float, 2U> &in1
//                const coder::array<float, 2U> &in2
//                const coder::array<float, 2U> &in3
// Return Type  : void
//
} // namespace checkerboard
} // namespace calibration
} // namespace internal
} // namespace vision
} // namespace coder
static void d_binary_expand_op(coder::array<float, 2U> &in1,
                               const coder::array<float, 2U> &in2,
                               const coder::array<float, 2U> &in3)
{
  coder::array<float, 2U> r;
  int aux_0_1;
  int aux_1_1;
  int aux_2_1;
  int c_loop_ub;
  int i;
  int i1;
  int i4;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  int stride_2_0;
  int stride_2_1;
  if (in3.size(0) == 1) {
    i = in2.size(0);
  } else {
    i = in3.size(0);
  }
  if (i == 1) {
    i = in1.size(0);
  } else if (in3.size(0) == 1) {
    i = in2.size(0);
  } else {
    i = in3.size(0);
  }
  if (in3.size(1) == 1) {
    i1 = in2.size(1);
  } else {
    i1 = in3.size(1);
  }
  if (i1 == 1) {
    i1 = in1.size(1);
  } else if (in3.size(1) == 1) {
    i1 = in2.size(1);
  } else {
    i1 = in3.size(1);
  }
  r.set_size(i, i1);
  stride_0_0 = (in1.size(0) != 1);
  stride_0_1 = (in1.size(1) != 1);
  stride_1_0 = (in2.size(0) != 1);
  stride_1_1 = (in2.size(1) != 1);
  stride_2_0 = (in3.size(0) != 1);
  stride_2_1 = (in3.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  aux_2_1 = 0;
  if (in3.size(1) == 1) {
    i = in2.size(1);
  } else {
    i = in3.size(1);
  }
  if (i == 1) {
    loop_ub = in1.size(1);
  } else if (in3.size(1) == 1) {
    loop_ub = in2.size(1);
  } else {
    loop_ub = in3.size(1);
  }
  for (i = 0; i < loop_ub; i++) {
    int b_loop_ub;
    int i3;
    i1 = in3.size(0);
    b_loop_ub = in2.size(0);
    if (i1 == 1) {
      i3 = b_loop_ub;
    } else {
      i3 = i1;
    }
    if (i3 == 1) {
      b_loop_ub = in1.size(0);
    } else if (i1 != 1) {
      b_loop_ub = i1;
    }
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      r[i1 + r.size(0) * i] =
          16.0F * in1[i1 * stride_0_0 + in1.size(0) * aux_0_1] -
          6.0F * (in2[i1 * stride_1_0 + in2.size(0) * aux_1_1] +
                  in3[i1 * stride_2_0 + in3.size(0) * aux_2_1]);
    }
    aux_2_1 += stride_2_1;
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  in1.set_size(r.size(0), r.size(1));
  loop_ub = r.size(1);
#pragma omp parallel for num_threads(32 > omp_get_max_threads()                \
                                         ? omp_get_max_threads()               \
                                         : 32) private(i4, c_loop_ub)

  for (int i2 = 0; i2 < loop_ub; i2++) {
    c_loop_ub = r.size(0);
    for (i4 = 0; i4 < c_loop_ub; i4++) {
      in1[i4 + in1.size(0) * i2] = r[i4 + r.size(0) * i2];
    }
  }
}

//
// Arguments    : const ::coder::array<float, 2U> &b_I
//                boolean_T highDistortion
//                ::coder::array<float, 2U> &cxy
//                ::coder::array<float, 2U> &c45
//                ::coder::array<float, 2U> &Ix
//                ::coder::array<float, 2U> &Iy
//                ::coder::array<float, 2U> &Ixy
//                ::coder::array<float, 2U> &I_45_45
// Return Type  : void
//
namespace coder {
namespace vision {
namespace internal {
namespace calibration {
namespace checkerboard {
void b_secondDerivCornerMetric(
    const ::coder::array<float, 2U> &b_I, boolean_T highDistortion,
    ::coder::array<float, 2U> &cxy, ::coder::array<float, 2U> &c45,
    ::coder::array<float, 2U> &Ix, ::coder::array<float, 2U> &Iy,
    ::coder::array<float, 2U> &Ixy, ::coder::array<float, 2U> &I_45_45)
{
  array<float, 3U> IfilterDdot;
  array<float, 3U> IfilterDot;
  array<float, 3U> b_IfilterDdot;
  array<float, 3U> x;
  array<float, 2U> Ig;
  array<float, 2U> checkerEdge;
  array<float, 2U> checkerEdgeComp;
  array<float, 2U> checkerboardPeaks;
  array<float, 2U> y;
  array<int, 1U> r;
  array<int, 1U> r1;
  int loop_ub;
  Ig.set_size(b_I.size(0), b_I.size(1));
  loop_ub = b_I.size(0) * b_I.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int k{0}; k < loop_ub; k++) {
      Ig[k] = b_I[k];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k < loop_ub; k++) {
      Ig[k] = b_I[k];
    }
  }
  d_imfilter(Ig);
  Iy.set_size(Ig.size(0), Ig.size(1));
  loop_ub = Ig.size(0) * Ig.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int k{0}; k < loop_ub; k++) {
      Iy[k] = Ig[k];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k < loop_ub; k++) {
      Iy[k] = Ig[k];
    }
  }
  imfilter(Iy);
  Ix.set_size(Ig.size(0), Ig.size(1));
  loop_ub = Ig.size(0) * Ig.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int k{0}; k < loop_ub; k++) {
      Ix[k] = Ig[k];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k < loop_ub; k++) {
      Ix[k] = Ig[k];
    }
  }
  b_imfilter(Ix);
  Ixy.set_size(Ix.size(0), Ix.size(1));
  loop_ub = Ix.size(0) * Ix.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int k{0}; k < loop_ub; k++) {
      Ixy[k] = Ix[k];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k < loop_ub; k++) {
      Ixy[k] = Ix[k];
    }
  }
  imfilter(Ixy);
  c45.set_size(b_I.size(0), b_I.size(1));
  loop_ub = b_I.size(0) * b_I.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int k{0}; k < loop_ub; k++) {
      c45[k] = 0.0F;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k < loop_ub; k++) {
      c45[k] = 0.0F;
    }
  }
  I_45_45.set_size(b_I.size(0), b_I.size(1));
  loop_ub = b_I.size(0) * b_I.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int k{0}; k < loop_ub; k++) {
      I_45_45[k] = 0.0F;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k < loop_ub; k++) {
      I_45_45[k] = 0.0F;
    }
  }
  if (!highDistortion) {
    int b_i;
    int e_loop_ub;
    int i1;
    int nx;
    int vstride;
    Ig.set_size(Ix.size(0), Ix.size(1));
    loop_ub = Ix.size(0) * Ix.size(1);
    if (static_cast<int>(loop_ub < 3200)) {
      for (int k{0}; k < loop_ub; k++) {
        Ig[k] = Ix[k] * 0.707106769F;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int k = 0; k < loop_ub; k++) {
        Ig[k] = Ix[k] * 0.707106769F;
      }
    }
    if ((Ig.size(0) == Iy.size(0)) && (Ig.size(1) == Iy.size(1))) {
      checkerEdgeComp.set_size(Ig.size(0), Ig.size(1));
      loop_ub = Ig.size(0) * Ig.size(1);
      if (static_cast<int>(loop_ub < 3200)) {
        for (int k{0}; k < loop_ub; k++) {
          checkerEdgeComp[k] = Ig[k] + Iy[k] * 0.707106769F;
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int k = 0; k < loop_ub; k++) {
          checkerEdgeComp[k] = Ig[k] + Iy[k] * 0.707106769F;
        }
      }
    } else {
      binary_expand_op(checkerEdgeComp, Ig, Iy);
    }
    I_45_45.set_size(checkerEdgeComp.size(0), checkerEdgeComp.size(1));
    loop_ub = checkerEdgeComp.size(0) * checkerEdgeComp.size(1);
    if (static_cast<int>(loop_ub < 3200)) {
      for (int k{0}; k < loop_ub; k++) {
        I_45_45[k] = checkerEdgeComp[k];
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int k = 0; k < loop_ub; k++) {
        I_45_45[k] = checkerEdgeComp[k];
      }
    }
    b_imfilter(I_45_45);
    checkerEdge.set_size(checkerEdgeComp.size(0), checkerEdgeComp.size(1));
    loop_ub = checkerEdgeComp.size(0) * checkerEdgeComp.size(1);
    if (static_cast<int>(loop_ub < 3200)) {
      for (int k{0}; k < loop_ub; k++) {
        checkerEdge[k] = checkerEdgeComp[k];
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int k = 0; k < loop_ub; k++) {
        checkerEdge[k] = checkerEdgeComp[k];
      }
    }
    imfilter(checkerEdge);
    if ((I_45_45.size(0) == checkerEdge.size(0)) &&
        (I_45_45.size(1) == checkerEdge.size(1))) {
      loop_ub = I_45_45.size(0) * I_45_45.size(1);
      if (static_cast<int>(loop_ub < 3200)) {
        for (int k{0}; k < loop_ub; k++) {
          I_45_45[k] =
              I_45_45[k] * 0.707106769F + checkerEdge[k] * -0.707106769F;
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int k = 0; k < loop_ub; k++) {
          I_45_45[k] =
              I_45_45[k] * 0.707106769F + checkerEdge[k] * -0.707106769F;
        }
      }
    } else {
      b_binary_expand_op(I_45_45, checkerEdge);
    }
    nx = Ixy.size(0) * Ixy.size(1);
    cxy.set_size(Ixy.size(0), Ixy.size(1));
    if (static_cast<int>(nx < 3200)) {
      for (int k{0}; k < nx; k++) {
        cxy[k] = std::abs(Ixy[k]);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int k = 0; k < nx; k++) {
        cxy[k] = std::abs(Ixy[k]);
      }
    }
    nx = checkerEdgeComp.size(0) * checkerEdgeComp.size(1);
    y.set_size(checkerEdgeComp.size(0), checkerEdgeComp.size(1));
    if (static_cast<int>(nx < 3200)) {
      for (int k{0}; k < nx; k++) {
        y[k] = std::abs(checkerEdgeComp[k]);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int k = 0; k < nx; k++) {
        y[k] = std::abs(checkerEdgeComp[k]);
      }
    }
    if ((Ig.size(0) == Iy.size(0)) && (Ig.size(1) == Iy.size(1))) {
      loop_ub = Ig.size(0) * Ig.size(1);
      if (static_cast<int>(loop_ub < 3200)) {
        for (int k{0}; k < loop_ub; k++) {
          Ig[k] = Ig[k] + Iy[k] * -0.707106769F;
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int k = 0; k < loop_ub; k++) {
          Ig[k] = Ig[k] + Iy[k] * -0.707106769F;
        }
      }
    } else {
      binary_expand_op(Ig, Iy);
    }
    nx = Ig.size(0) * Ig.size(1);
    checkerEdge.set_size(Ig.size(0), Ig.size(1));
    if (static_cast<int>(nx < 3200)) {
      for (int k{0}; k < nx; k++) {
        checkerEdge[k] = std::abs(Ig[k]);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int k = 0; k < nx; k++) {
        checkerEdge[k] = std::abs(Ig[k]);
      }
    }
    if (y.size(0) == 1) {
      b_i = checkerEdge.size(0);
    } else {
      b_i = y.size(0);
    }
    if (y.size(1) == 1) {
      i1 = checkerEdge.size(1);
    } else {
      i1 = y.size(1);
    }
    if ((y.size(0) == checkerEdge.size(0)) &&
        (y.size(1) == checkerEdge.size(1)) && (cxy.size(0) == b_i) &&
        (cxy.size(1) == i1)) {
      loop_ub = cxy.size(0) * cxy.size(1);
      if (static_cast<int>(loop_ub < 3200)) {
        for (int k{0}; k < loop_ub; k++) {
          cxy[k] = 16.0F * cxy[k] - 6.0F * (y[k] + checkerEdge[k]);
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int k = 0; k < loop_ub; k++) {
          cxy[k] = 16.0F * cxy[k] - 6.0F * (y[k] + checkerEdge[k]);
        }
      }
    } else {
      d_binary_expand_op(cxy, y, checkerEdge);
    }
    e_loop_ub = cxy.size(0) * cxy.size(1) - 1;
    vstride = 0;
    for (int i{0}; i <= e_loop_ub; i++) {
      if (cxy[i] < 0.0F) {
        vstride++;
      }
    }
    r.set_size(vstride);
    vstride = 0;
    for (int i{0}; i <= e_loop_ub; i++) {
      if (cxy[i] < 0.0F) {
        r[vstride] = i + 1;
        vstride++;
      }
    }
    loop_ub = r.size(0);
    for (int i{0}; i < loop_ub; i++) {
      cxy[r[i] - 1] = 0.0F;
    }
    nx = I_45_45.size(0) * I_45_45.size(1);
    c45.set_size(I_45_45.size(0), I_45_45.size(1));
    if (static_cast<int>(nx < 3200)) {
      for (int k{0}; k < nx; k++) {
        c45[k] = std::abs(I_45_45[k]);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int k = 0; k < nx; k++) {
        c45[k] = std::abs(I_45_45[k]);
      }
    }
    nx = Ix.size(0) * Ix.size(1);
    y.set_size(Ix.size(0), Ix.size(1));
    if (static_cast<int>(nx < 3200)) {
      for (int k{0}; k < nx; k++) {
        y[k] = std::abs(Ix[k]);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int k = 0; k < nx; k++) {
        y[k] = std::abs(Ix[k]);
      }
    }
    nx = Iy.size(0) * Iy.size(1);
    checkerEdge.set_size(Iy.size(0), Iy.size(1));
    if (static_cast<int>(nx < 3200)) {
      for (int k{0}; k < nx; k++) {
        checkerEdge[k] = std::abs(Iy[k]);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int k = 0; k < nx; k++) {
        checkerEdge[k] = std::abs(Iy[k]);
      }
    }
    if (y.size(0) == 1) {
      b_i = checkerEdge.size(0);
    } else {
      b_i = y.size(0);
    }
    if (y.size(1) == 1) {
      i1 = checkerEdge.size(1);
    } else {
      i1 = y.size(1);
    }
    if ((y.size(0) == checkerEdge.size(0)) &&
        (y.size(1) == checkerEdge.size(1)) && (c45.size(0) == b_i) &&
        (c45.size(1) == i1)) {
      loop_ub = c45.size(0) * c45.size(1);
      if (static_cast<int>(loop_ub < 3200)) {
        for (int k{0}; k < loop_ub; k++) {
          c45[k] = 16.0F * c45[k] - 6.0F * (y[k] + checkerEdge[k]);
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int k = 0; k < loop_ub; k++) {
          c45[k] = 16.0F * c45[k] - 6.0F * (y[k] + checkerEdge[k]);
        }
      }
    } else {
      d_binary_expand_op(c45, y, checkerEdge);
    }
    e_loop_ub = c45.size(0) * c45.size(1) - 1;
    vstride = 0;
    for (int i{0}; i <= e_loop_ub; i++) {
      if (c45[i] < 0.0F) {
        vstride++;
      }
    }
    r1.set_size(vstride);
    vstride = 0;
    for (int i{0}; i <= e_loop_ub; i++) {
      if (c45[i] < 0.0F) {
        r1[vstride] = i + 1;
        vstride++;
      }
    }
    loop_ub = r1.size(0);
    for (int i{0}; i < loop_ub; i++) {
      c45[r1[i] - 1] = 0.0F;
    }
  } else {
    cell_wrap_11 derivFilters[12];
    int b_loop_ub;
    int b_nx;
    int c_loop_ub;
    int d_loop_ub;
    int nx;
    getDerivFilters(derivFilters);
    checkerboardPeaks.set_size(b_I.size(0), b_I.size(1));
    loop_ub = b_I.size(0) * b_I.size(1);
    if (static_cast<int>(loop_ub < 3200)) {
      for (int k{0}; k < loop_ub; k++) {
        checkerboardPeaks[k] = 0.0F;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int k = 0; k < loop_ub; k++) {
        checkerboardPeaks[k] = 0.0F;
      }
    }
    IfilterDot.set_size(Ig.size(0), Ig.size(1), 4);
    loop_ub = (Ig.size(0) * Ig.size(1)) << 2;
    if (static_cast<int>(loop_ub < 3200)) {
      for (int k{0}; k < loop_ub; k++) {
        IfilterDot[k] = 0.0F;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int k = 0; k < loop_ub; k++) {
        IfilterDot[k] = 0.0F;
      }
    }
    IfilterDdot.set_size(Ig.size(0), Ig.size(1), 4);
    loop_ub = (Ig.size(0) * Ig.size(1)) << 2;
    if (static_cast<int>(loop_ub < 3200)) {
      for (int k{0}; k < loop_ub; k++) {
        IfilterDdot[k] = 0.0F;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int k = 0; k < loop_ub; k++) {
        IfilterDdot[k] = 0.0F;
      }
    }
    loop_ub = Ig.size(0) * Ig.size(1);
    b_loop_ub = Ig.size(0) * Ig.size(1);
    c_loop_ub = Ig.size(0) * Ig.size(1);
    d_loop_ub = Ig.size(0) * Ig.size(1);
    nx = Ix.size(0) * Ix.size(1);
    b_nx = Iy.size(0) * Iy.size(1);
    for (int filtIdx{0}; filtIdx < 3; filtIdx++) {
      float varargin_1;
      int b_i;
      int e_loop_ub;
      int i1;
      int vstride;
      checkerEdge.set_size(Ig.size(0), Ig.size(1));
      for (int i{0}; i < loop_ub; i++) {
        checkerEdge[i] = Ig[i];
      }
      imfilter(checkerEdge, derivFilters[filtIdx].f1);
      vstride = checkerEdge.size(1);
      for (int i{0}; i < vstride; i++) {
        e_loop_ub = checkerEdge.size(0);
        for (b_i = 0; b_i < e_loop_ub; b_i++) {
          IfilterDot[b_i + IfilterDot.size(0) * i] =
              checkerEdge[b_i + checkerEdge.size(0) * i];
        }
      }
      checkerEdge.set_size(Ig.size(0), Ig.size(1));
      for (int i{0}; i < b_loop_ub; i++) {
        checkerEdge[i] = Ig[i];
      }
      imfilter(checkerEdge, derivFilters[filtIdx + 3].f1);
      vstride = checkerEdge.size(1);
      for (int i{0}; i < vstride; i++) {
        e_loop_ub = checkerEdge.size(0);
        for (b_i = 0; b_i < e_loop_ub; b_i++) {
          IfilterDot[(b_i + IfilterDot.size(0) * i) +
                     IfilterDot.size(0) * IfilterDot.size(1)] =
              checkerEdge[b_i + checkerEdge.size(0) * i];
        }
      }
      checkerEdge.set_size(Ig.size(0), Ig.size(1));
      for (int i{0}; i < c_loop_ub; i++) {
        checkerEdge[i] = Ig[i];
      }
      imfilter(checkerEdge, derivFilters[filtIdx + 6].f1);
      vstride = checkerEdge.size(1);
      for (int i{0}; i < vstride; i++) {
        e_loop_ub = checkerEdge.size(0);
        for (b_i = 0; b_i < e_loop_ub; b_i++) {
          IfilterDot[(b_i + IfilterDot.size(0) * i) +
                     IfilterDot.size(0) * IfilterDot.size(1) * 2] =
              checkerEdge[b_i + checkerEdge.size(0) * i];
        }
      }
      checkerEdge.set_size(Ig.size(0), Ig.size(1));
      for (int i{0}; i < d_loop_ub; i++) {
        checkerEdge[i] = Ig[i];
      }
      imfilter(checkerEdge, derivFilters[filtIdx + 9].f1);
      vstride = checkerEdge.size(1);
      for (int i{0}; i < vstride; i++) {
        e_loop_ub = checkerEdge.size(0);
        for (b_i = 0; b_i < e_loop_ub; b_i++) {
          IfilterDot[(b_i + IfilterDot.size(0) * i) +
                     IfilterDot.size(0) * IfilterDot.size(1) * 3] =
              checkerEdge[b_i + checkerEdge.size(0) * i];
        }
      }
      vstride = IfilterDot.size(1);
      checkerEdge.set_size(IfilterDot.size(0), IfilterDot.size(1));
      for (int i{0}; i < vstride; i++) {
        e_loop_ub = IfilterDot.size(0);
        for (b_i = 0; b_i < e_loop_ub; b_i++) {
          checkerEdge[b_i + checkerEdge.size(0) * i] =
              IfilterDot[b_i + IfilterDot.size(0) * i];
        }
      }
      imfilter(checkerEdge, derivFilters[filtIdx].f1);
      vstride = checkerEdge.size(1);
      for (int i{0}; i < vstride; i++) {
        e_loop_ub = checkerEdge.size(0);
        for (b_i = 0; b_i < e_loop_ub; b_i++) {
          IfilterDdot[b_i + IfilterDdot.size(0) * i] =
              checkerEdge[b_i + checkerEdge.size(0) * i];
        }
      }
      vstride = IfilterDot.size(1);
      checkerEdge.set_size(IfilterDot.size(0), IfilterDot.size(1));
      for (int i{0}; i < vstride; i++) {
        e_loop_ub = IfilterDot.size(0);
        for (b_i = 0; b_i < e_loop_ub; b_i++) {
          checkerEdge[b_i + checkerEdge.size(0) * i] =
              IfilterDot[(b_i + IfilterDot.size(0) * i) +
                         IfilterDot.size(0) * IfilterDot.size(1)];
        }
      }
      imfilter(checkerEdge, derivFilters[filtIdx + 3].f1);
      vstride = checkerEdge.size(1);
      for (int i{0}; i < vstride; i++) {
        e_loop_ub = checkerEdge.size(0);
        for (b_i = 0; b_i < e_loop_ub; b_i++) {
          IfilterDdot[(b_i + IfilterDdot.size(0) * i) +
                      IfilterDdot.size(0) * IfilterDdot.size(1)] =
              checkerEdge[b_i + checkerEdge.size(0) * i];
        }
      }
      vstride = IfilterDot.size(1);
      checkerEdge.set_size(IfilterDot.size(0), IfilterDot.size(1));
      for (int i{0}; i < vstride; i++) {
        e_loop_ub = IfilterDot.size(0);
        for (b_i = 0; b_i < e_loop_ub; b_i++) {
          checkerEdge[b_i + checkerEdge.size(0) * i] =
              IfilterDot[b_i + IfilterDot.size(0) * i];
        }
      }
      imfilter(checkerEdge, derivFilters[filtIdx + 6].f1);
      vstride = checkerEdge.size(1);
      for (int i{0}; i < vstride; i++) {
        e_loop_ub = checkerEdge.size(0);
        for (b_i = 0; b_i < e_loop_ub; b_i++) {
          IfilterDdot[(b_i + IfilterDdot.size(0) * i) +
                      IfilterDdot.size(0) * IfilterDdot.size(1) * 2] =
              checkerEdge[b_i + checkerEdge.size(0) * i];
        }
      }
      vstride = IfilterDot.size(1);
      checkerEdge.set_size(IfilterDot.size(0), IfilterDot.size(1));
      for (int i{0}; i < vstride; i++) {
        e_loop_ub = IfilterDot.size(0);
        for (b_i = 0; b_i < e_loop_ub; b_i++) {
          checkerEdge[b_i + checkerEdge.size(0) * i] =
              IfilterDot[(b_i + IfilterDot.size(0) * i) +
                         IfilterDot.size(0) * IfilterDot.size(1)];
        }
      }
      imfilter(checkerEdge, derivFilters[filtIdx + 9].f1);
      vstride = checkerEdge.size(1);
      for (int i{0}; i < vstride; i++) {
        e_loop_ub = checkerEdge.size(0);
        for (b_i = 0; b_i < e_loop_ub; b_i++) {
          IfilterDdot[(b_i + IfilterDdot.size(0) * i) +
                      IfilterDdot.size(0) * IfilterDdot.size(1) * 3] =
              checkerEdge[b_i + checkerEdge.size(0) * i];
        }
      }
      checkerEdgeComp.set_size(Ix.size(0), Ix.size(1));
      for (vstride = 0; vstride < nx; vstride++) {
        checkerEdgeComp[vstride] = std::abs(Ix[vstride]);
      }
      y.set_size(Iy.size(0), Iy.size(1));
      for (vstride = 0; vstride < b_nx; vstride++) {
        y[vstride] = std::abs(Iy[vstride]);
      }
      if ((checkerEdgeComp.size(0) == y.size(0)) &&
          (checkerEdgeComp.size(1) == y.size(1))) {
        vstride = checkerEdgeComp.size(0) * checkerEdgeComp.size(1);
        for (int i{0}; i < vstride; i++) {
          checkerEdgeComp[i] = checkerEdgeComp[i] * 0.3F + y[i] * 0.3F;
        }
      } else {
        c_binary_expand_op(checkerEdgeComp, y);
      }
      vstride = IfilterDdot.size(1);
      if ((IfilterDdot.size(0) == checkerEdgeComp.size(0)) &&
          (IfilterDdot.size(1) == checkerEdgeComp.size(1))) {
        checkerEdge.set_size(IfilterDdot.size(0), IfilterDdot.size(1));
        for (int i{0}; i < vstride; i++) {
          e_loop_ub = IfilterDdot.size(0);
          for (b_i = 0; b_i < e_loop_ub; b_i++) {
            checkerEdge[b_i + checkerEdge.size(0) * i] =
                (IfilterDdot[b_i + IfilterDdot.size(0) * i] +
                 IfilterDdot[(b_i + IfilterDdot.size(0) * i) +
                             IfilterDdot.size(0) * IfilterDdot.size(1)]) -
                checkerEdgeComp[b_i + checkerEdgeComp.size(0) * i];
          }
        }
      } else {
        binary_expand_op(checkerEdge, IfilterDdot, checkerEdgeComp);
      }
      vstride = IfilterDdot.size(1);
      if ((IfilterDdot.size(0) == checkerEdgeComp.size(0)) &&
          (IfilterDdot.size(1) == checkerEdgeComp.size(1))) {
        checkerEdgeComp.set_size(IfilterDdot.size(0), IfilterDdot.size(1));
        for (int i{0}; i < vstride; i++) {
          e_loop_ub = IfilterDdot.size(0);
          for (b_i = 0; b_i < e_loop_ub; b_i++) {
            checkerEdgeComp[b_i + checkerEdgeComp.size(0) * i] =
                (IfilterDdot[(b_i + IfilterDdot.size(0) * i) +
                             IfilterDdot.size(0) * IfilterDdot.size(1) * 2] +
                 IfilterDdot[(b_i + IfilterDdot.size(0) * i) +
                             IfilterDdot.size(0) * IfilterDdot.size(1) * 3]) -
                checkerEdgeComp[b_i + checkerEdgeComp.size(0) * i];
          }
        }
      } else {
        binary_expand_op(checkerEdgeComp, IfilterDdot);
      }
      e_loop_ub = checkerEdge.size(0) * checkerEdge.size(1) - 1;
      vstride = 0;
      for (int i{0}; i <= e_loop_ub; i++) {
        if (checkerEdge[i] < 0.0F) {
          vstride++;
        }
      }
      r.set_size(vstride);
      vstride = 0;
      for (int i{0}; i <= e_loop_ub; i++) {
        if (checkerEdge[i] < 0.0F) {
          r[vstride] = i + 1;
          vstride++;
        }
      }
      vstride = r.size(0);
      for (int i{0}; i < vstride; i++) {
        checkerEdge[r[i] - 1] = 0.0F;
      }
      e_loop_ub = checkerEdgeComp.size(0) * checkerEdgeComp.size(1) - 1;
      vstride = 0;
      for (int i{0}; i <= e_loop_ub; i++) {
        if (checkerEdgeComp[i] < 0.0F) {
          vstride++;
        }
      }
      r1.set_size(vstride);
      vstride = 0;
      for (int i{0}; i <= e_loop_ub; i++) {
        if (checkerEdgeComp[i] < 0.0F) {
          r1[vstride] = i + 1;
          vstride++;
        }
      }
      vstride = r1.size(0);
      for (int i{0}; i < vstride; i++) {
        checkerEdgeComp[r1[i] - 1] = 0.0F;
      }
      y.set_size(checkerEdge.size(0), checkerEdge.size(1));
      vstride = checkerEdge.size(0) * checkerEdge.size(1) - 1;
      for (int i{0}; i <= vstride; i++) {
        y[i] = checkerEdge[i];
      }
      imdilate(y, checkerEdge);
      y.set_size(checkerEdgeComp.size(0), checkerEdgeComp.size(1));
      vstride = checkerEdgeComp.size(0) * checkerEdgeComp.size(1) - 1;
      for (int i{0}; i <= vstride; i++) {
        y[i] = checkerEdgeComp[i];
      }
      imdilate(y, checkerEdgeComp);
      vstride = IfilterDdot.size(1);
      b_IfilterDdot.set_size(IfilterDdot.size(0), IfilterDdot.size(1), 2);
      for (int i{0}; i < 2; i++) {
        for (b_i = 0; b_i < vstride; b_i++) {
          e_loop_ub = IfilterDdot.size(0);
          for (i1 = 0; i1 < e_loop_ub; i1++) {
            b_IfilterDdot[(i1 + b_IfilterDdot.size(0) * b_i) +
                          b_IfilterDdot.size(0) * b_IfilterDdot.size(1) * i] =
                IfilterDdot[(i1 + IfilterDdot.size(0) * b_i) +
                            IfilterDdot.size(0) * IfilterDdot.size(1) * i];
          }
        }
      }
      bsxfun(b_IfilterDdot, checkerEdge, x);
      vstride = x.size(0) * x.size(1) * 2;
      x.set_size(x.size(0), x.size(1), 2);
      for (int i{0}; i < vstride; i++) {
        varargin_1 = x[i];
        x[i] = std::fmax(varargin_1, 0.0F);
      }
      if ((x.size(0) == 0) || (x.size(1) == 0)) {
        y.set_size(x.size(0), x.size(1));
        vstride = x.size(0) * x.size(1);
        for (int i{0}; i < vstride; i++) {
          y[i] = 0.0F;
        }
      } else {
        vstride = x.size(0) * x.size(1);
        y.set_size(x.size(0), x.size(1));
        for (e_loop_ub = 0; e_loop_ub < vstride; e_loop_ub++) {
          y[e_loop_ub] = x[e_loop_ub];
        }
        for (e_loop_ub = 0; e_loop_ub < vstride; e_loop_ub++) {
          y[e_loop_ub] = y[e_loop_ub] + x[vstride + e_loop_ub];
        }
      }
      vstride = IfilterDdot.size(1);
      b_IfilterDdot.set_size(IfilterDdot.size(0), IfilterDdot.size(1), 2);
      for (int i{0}; i < 2; i++) {
        for (b_i = 0; b_i < vstride; b_i++) {
          e_loop_ub = IfilterDdot.size(0);
          for (i1 = 0; i1 < e_loop_ub; i1++) {
            b_IfilterDdot[(i1 + b_IfilterDdot.size(0) * b_i) +
                          b_IfilterDdot.size(0) * b_IfilterDdot.size(1) * i] =
                IfilterDdot[(i1 + IfilterDdot.size(0) * b_i) +
                            IfilterDdot.size(0) * IfilterDdot.size(1) *
                                (i + 2)];
          }
        }
      }
      bsxfun(b_IfilterDdot, checkerEdgeComp, x);
      vstride = x.size(0) * x.size(1) * 2;
      x.set_size(x.size(0), x.size(1), 2);
      for (int i{0}; i < vstride; i++) {
        varargin_1 = x[i];
        x[i] = std::fmax(varargin_1, 0.0F);
      }
      if ((x.size(0) == 0) || (x.size(1) == 0)) {
        checkerEdge.set_size(x.size(0), x.size(1));
        vstride = x.size(0) * x.size(1);
        for (int i{0}; i < vstride; i++) {
          checkerEdge[i] = 0.0F;
        }
      } else {
        vstride = x.size(0) * x.size(1);
        checkerEdge.set_size(x.size(0), x.size(1));
        for (e_loop_ub = 0; e_loop_ub < vstride; e_loop_ub++) {
          checkerEdge[e_loop_ub] = x[e_loop_ub];
        }
        for (e_loop_ub = 0; e_loop_ub < vstride; e_loop_ub++) {
          checkerEdge[e_loop_ub] =
              checkerEdge[e_loop_ub] + x[vstride + e_loop_ub];
        }
      }
      if (checkerboardPeaks.size(0) == 1) {
        b_i = y.size(0);
      } else {
        b_i = checkerboardPeaks.size(0);
      }
      if (checkerboardPeaks.size(1) == 1) {
        i1 = y.size(1);
      } else {
        i1 = checkerboardPeaks.size(1);
      }
      if ((checkerboardPeaks.size(0) == y.size(0)) &&
          (checkerboardPeaks.size(1) == y.size(1)) &&
          (b_i == checkerEdge.size(0)) && (i1 == checkerEdge.size(1))) {
        vstride = checkerboardPeaks.size(0) * checkerboardPeaks.size(1);
        for (int i{0}; i < vstride; i++) {
          checkerboardPeaks[i] = (checkerboardPeaks[i] + y[i]) + checkerEdge[i];
        }
      } else {
        b_binary_expand_op(checkerboardPeaks, y, checkerEdge);
      }
    }
    imerode(checkerboardPeaks, cxy);
  }
}

//
// Arguments    : const ::coder::array<float, 2U> &b_I
//                double sigma
//                boolean_T highDistortion
//                ::coder::array<float, 2U> &cxy
//                ::coder::array<float, 2U> &c45
//                ::coder::array<float, 2U> &Ix
//                ::coder::array<float, 2U> &Iy
//                ::coder::array<float, 2U> &Ixy
//                ::coder::array<float, 2U> &I_45_45
// Return Type  : void
//
void secondDerivCornerMetric(const ::coder::array<float, 2U> &b_I, double sigma,
                             boolean_T highDistortion,
                             ::coder::array<float, 2U> &cxy,
                             ::coder::array<float, 2U> &c45,
                             ::coder::array<float, 2U> &Ix,
                             ::coder::array<float, 2U> &Iy,
                             ::coder::array<float, 2U> &Ixy,
                             ::coder::array<float, 2U> &I_45_45)
{
  array<double, 2U> b_y;
  array<double, 2U> y;
  array<double, 1U> b_G_data;
  array<double, 1U> c_G_data;
  array<float, 3U> IfilterDdot;
  array<float, 3U> IfilterDot;
  array<float, 3U> b_IfilterDdot;
  array<float, 3U> x;
  array<float, 2U> Ig;
  array<float, 2U> c_y;
  array<float, 2U> checkerEdge;
  array<float, 2U> checkerEdgeComp;
  array<float, 2U> checkerboardPeaks;
  array<int, 1U> r;
  array<int, 1U> r1;
  double G_data[225];
  double y_data[225];
  double a;
  double siz_idx_0;
  int G_size[2];
  int y_size[2];
  int b_i;
  int i;
  int k;
  int loop_ub;
  int nx;
  int ny;
  int partialTrueCount;
  int trueCount;
  unsigned char b_tmp_data[225];
  boolean_T tmp_data[225];
  siz_idx_0 = ((std::round(sigma * 7.0) + 1.0) - 1.0) / 2.0;
  if (std::isnan(-siz_idx_0) || std::isnan(siz_idx_0)) {
    y.set_size(1, 1);
    y[0] = rtNaN;
  } else if (siz_idx_0 < -siz_idx_0) {
    y.set_size(1, 0);
  } else if ((std::isinf(-siz_idx_0) || std::isinf(siz_idx_0)) &&
             (-siz_idx_0 == siz_idx_0)) {
    y.set_size(1, 1);
    y[0] = rtNaN;
  } else if (std::floor(-siz_idx_0) == -siz_idx_0) {
    a = -siz_idx_0;
    loop_ub = static_cast<int>(siz_idx_0 - (-siz_idx_0));
    y.set_size(1, loop_ub + 1);
    if (static_cast<int>(loop_ub + 1 < 3200)) {
      for (int j{0}; j <= loop_ub; j++) {
        y[j] = -siz_idx_0 + static_cast<double>(j);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int j = 0; j <= loop_ub; j++) {
        y[j] = a + static_cast<double>(j);
      }
    }
  } else {
    eml_float_colon(-siz_idx_0, siz_idx_0, y);
  }
  if (std::isnan(-siz_idx_0) || std::isnan(siz_idx_0)) {
    b_y.set_size(1, 1);
    b_y[0] = rtNaN;
  } else if (siz_idx_0 < -siz_idx_0) {
    b_y.set_size(1, 0);
  } else if ((std::isinf(-siz_idx_0) || std::isinf(siz_idx_0)) &&
             (-siz_idx_0 == siz_idx_0)) {
    b_y.set_size(1, 1);
    b_y[0] = rtNaN;
  } else if (std::floor(-siz_idx_0) == -siz_idx_0) {
    a = -siz_idx_0;
    loop_ub = static_cast<int>(siz_idx_0 - (-siz_idx_0));
    b_y.set_size(1, loop_ub + 1);
    if (static_cast<int>(loop_ub + 1 < 3200)) {
      for (int j{0}; j <= loop_ub; j++) {
        b_y[j] = -siz_idx_0 + static_cast<double>(j);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int j = 0; j <= loop_ub; j++) {
        b_y[j] = a + static_cast<double>(j);
      }
    }
  } else {
    eml_float_colon(-siz_idx_0, siz_idx_0, b_y);
  }
  nx = y.size(1);
  ny = b_y.size(1);
  G_size[0] = b_y.size(1);
  G_size[1] = y.size(1);
  y_size[0] = b_y.size(1);
  y_size[1] = y.size(1);
  if (static_cast<int>(y.size(1) * b_y.size(1) < 3200)) {
    for (int j{0}; j < nx; j++) {
      for (i = 0; i < ny; i++) {
        G_data[i + G_size[0] * j] = y[j];
        y_data[i + y_size[0] * j] = b_y[i];
      }
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32) private(i)

    for (int j = 0; j < nx; j++) {
      for (i = 0; i < ny; i++) {
        G_data[i + G_size[0] * j] = y[j];
        y_data[i + y_size[0] * j] = b_y[i];
      }
    }
  }
  if ((G_size[0] == y_size[0]) && (G_size[1] == y_size[1])) {
    a = 2.0 * sigma * sigma;
    loop_ub = G_size[0] * G_size[1];
    for (b_i = 0; b_i < loop_ub; b_i++) {
      siz_idx_0 = y_data[b_i];
      G_data[b_i] = -(G_data[b_i] * G_data[b_i] + siz_idx_0 * siz_idx_0) / a;
    }
  } else {
    binary_expand_op(G_data, G_size, y_data, y_size, sigma);
  }
  ny = G_size[0] * G_size[1];
  for (k = 0; k < ny; k++) {
    G_data[k] = std::exp(G_data[k]);
  }
  b_G_data.set(&G_data[0], ny);
  a = 2.2204460492503131E-16 * ::coder::internal::maximum(b_G_data);
  for (b_i = 0; b_i < ny; b_i++) {
    tmp_data[b_i] = (G_data[b_i] < a);
  }
  k = G_size[0] * G_size[1] - 1;
  trueCount = 0;
  partialTrueCount = 0;
  for (b_i = 0; b_i <= k; b_i++) {
    if (tmp_data[b_i]) {
      trueCount++;
      b_tmp_data[partialTrueCount] = static_cast<unsigned char>(b_i + 1);
      partialTrueCount++;
    }
  }
  for (b_i = 0; b_i < trueCount; b_i++) {
    G_data[b_tmp_data[b_i] - 1] = 0.0;
  }
  c_G_data.set(&G_data[0], ny);
  siz_idx_0 = combineVectorElements(c_G_data);
  if (siz_idx_0 != 0.0) {
    for (b_i = 0; b_i < ny; b_i++) {
      G_data[b_i] /= siz_idx_0;
    }
  }
  Ig.set_size(b_I.size(0), b_I.size(1));
  loop_ub = b_I.size(0) * b_I.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int j{0}; j < loop_ub; j++) {
      Ig[j] = b_I[j];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int j = 0; j < loop_ub; j++) {
      Ig[j] = b_I[j];
    }
  }
  imfilter(Ig, G_data, G_size);
  Iy.set_size(Ig.size(0), Ig.size(1));
  loop_ub = Ig.size(0) * Ig.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int j{0}; j < loop_ub; j++) {
      Iy[j] = Ig[j];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int j = 0; j < loop_ub; j++) {
      Iy[j] = Ig[j];
    }
  }
  imfilter(Iy);
  Ix.set_size(Ig.size(0), Ig.size(1));
  loop_ub = Ig.size(0) * Ig.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int j{0}; j < loop_ub; j++) {
      Ix[j] = Ig[j];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int j = 0; j < loop_ub; j++) {
      Ix[j] = Ig[j];
    }
  }
  b_imfilter(Ix);
  Ixy.set_size(Ix.size(0), Ix.size(1));
  loop_ub = Ix.size(0) * Ix.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int j{0}; j < loop_ub; j++) {
      Ixy[j] = Ix[j];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int j = 0; j < loop_ub; j++) {
      Ixy[j] = Ix[j];
    }
  }
  imfilter(Ixy);
  c45.set_size(b_I.size(0), b_I.size(1));
  loop_ub = b_I.size(0) * b_I.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int j{0}; j < loop_ub; j++) {
      c45[j] = 0.0F;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int j = 0; j < loop_ub; j++) {
      c45[j] = 0.0F;
    }
  }
  I_45_45.set_size(b_I.size(0), b_I.size(1));
  loop_ub = b_I.size(0) * b_I.size(1);
  if (static_cast<int>(loop_ub < 3200)) {
    for (int j{0}; j < loop_ub; j++) {
      I_45_45[j] = 0.0F;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int j = 0; j < loop_ub; j++) {
      I_45_45[j] = 0.0F;
    }
  }
  if (!highDistortion) {
    Ig.set_size(Ix.size(0), Ix.size(1));
    loop_ub = Ix.size(0) * Ix.size(1);
    if (static_cast<int>(loop_ub < 3200)) {
      for (int j{0}; j < loop_ub; j++) {
        Ig[j] = Ix[j] * 0.707106769F;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int j = 0; j < loop_ub; j++) {
        Ig[j] = Ix[j] * 0.707106769F;
      }
    }
    if ((Ig.size(0) == Iy.size(0)) && (Ig.size(1) == Iy.size(1))) {
      checkerEdgeComp.set_size(Ig.size(0), Ig.size(1));
      loop_ub = Ig.size(0) * Ig.size(1);
      if (static_cast<int>(loop_ub < 3200)) {
        for (int j{0}; j < loop_ub; j++) {
          checkerEdgeComp[j] = Ig[j] + Iy[j] * 0.707106769F;
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int j = 0; j < loop_ub; j++) {
          checkerEdgeComp[j] = Ig[j] + Iy[j] * 0.707106769F;
        }
      }
    } else {
      binary_expand_op(checkerEdgeComp, Ig, Iy);
    }
    I_45_45.set_size(checkerEdgeComp.size(0), checkerEdgeComp.size(1));
    loop_ub = checkerEdgeComp.size(0) * checkerEdgeComp.size(1);
    if (static_cast<int>(loop_ub < 3200)) {
      for (int j{0}; j < loop_ub; j++) {
        I_45_45[j] = checkerEdgeComp[j];
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int j = 0; j < loop_ub; j++) {
        I_45_45[j] = checkerEdgeComp[j];
      }
    }
    b_imfilter(I_45_45);
    checkerEdge.set_size(checkerEdgeComp.size(0), checkerEdgeComp.size(1));
    loop_ub = checkerEdgeComp.size(0) * checkerEdgeComp.size(1);
    if (static_cast<int>(loop_ub < 3200)) {
      for (int j{0}; j < loop_ub; j++) {
        checkerEdge[j] = checkerEdgeComp[j];
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int j = 0; j < loop_ub; j++) {
        checkerEdge[j] = checkerEdgeComp[j];
      }
    }
    imfilter(checkerEdge);
    if ((I_45_45.size(0) == checkerEdge.size(0)) &&
        (I_45_45.size(1) == checkerEdge.size(1))) {
      loop_ub = I_45_45.size(0) * I_45_45.size(1);
      if (static_cast<int>(loop_ub < 3200)) {
        for (int j{0}; j < loop_ub; j++) {
          I_45_45[j] =
              I_45_45[j] * 0.707106769F + checkerEdge[j] * -0.707106769F;
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int j = 0; j < loop_ub; j++) {
          I_45_45[j] =
              I_45_45[j] * 0.707106769F + checkerEdge[j] * -0.707106769F;
        }
      }
    } else {
      b_binary_expand_op(I_45_45, checkerEdge);
    }
    siz_idx_0 = sigma * sigma;
    a = 1.5 * sigma;
    nx = Ixy.size(0) * Ixy.size(1);
    cxy.set_size(Ixy.size(0), Ixy.size(1));
    if (static_cast<int>(nx < 3200)) {
      for (int j{0}; j < nx; j++) {
        cxy[j] = std::abs(Ixy[j]);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int j = 0; j < nx; j++) {
        cxy[j] = std::abs(Ixy[j]);
      }
    }
    nx = checkerEdgeComp.size(0) * checkerEdgeComp.size(1);
    c_y.set_size(checkerEdgeComp.size(0), checkerEdgeComp.size(1));
    if (static_cast<int>(nx < 3200)) {
      for (int j{0}; j < nx; j++) {
        c_y[j] = std::abs(checkerEdgeComp[j]);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int j = 0; j < nx; j++) {
        c_y[j] = std::abs(checkerEdgeComp[j]);
      }
    }
    if ((Ig.size(0) == Iy.size(0)) && (Ig.size(1) == Iy.size(1))) {
      loop_ub = Ig.size(0) * Ig.size(1);
      if (static_cast<int>(loop_ub < 3200)) {
        for (int j{0}; j < loop_ub; j++) {
          Ig[j] = Ig[j] + Iy[j] * -0.707106769F;
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int j = 0; j < loop_ub; j++) {
          Ig[j] = Ig[j] + Iy[j] * -0.707106769F;
        }
      }
    } else {
      binary_expand_op(Ig, Iy);
    }
    nx = Ig.size(0) * Ig.size(1);
    checkerEdge.set_size(Ig.size(0), Ig.size(1));
    if (static_cast<int>(nx < 3200)) {
      for (int j{0}; j < nx; j++) {
        checkerEdge[j] = std::abs(Ig[j]);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int j = 0; j < nx; j++) {
        checkerEdge[j] = std::abs(Ig[j]);
      }
    }
    if (c_y.size(0) == 1) {
      ny = checkerEdge.size(0);
    } else {
      ny = c_y.size(0);
    }
    if (c_y.size(1) == 1) {
      b_i = checkerEdge.size(1);
    } else {
      b_i = c_y.size(1);
    }
    if ((c_y.size(0) == checkerEdge.size(0)) &&
        (c_y.size(1) == checkerEdge.size(1)) && (cxy.size(0) == ny) &&
        (cxy.size(1) == b_i)) {
      loop_ub = cxy.size(0) * cxy.size(1);
      if (static_cast<int>(loop_ub < 3200)) {
        for (int j{0}; j < loop_ub; j++) {
          cxy[j] = static_cast<float>(siz_idx_0) * cxy[j] -
                   static_cast<float>(a) * (c_y[j] + checkerEdge[j]);
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int j = 0; j < loop_ub; j++) {
          cxy[j] = static_cast<float>(siz_idx_0) * cxy[j] -
                   static_cast<float>(a) * (c_y[j] + checkerEdge[j]);
        }
      }
    } else {
      binary_expand_op(cxy, siz_idx_0, a, c_y, checkerEdge);
    }
    k = cxy.size(0) * cxy.size(1) - 1;
    trueCount = 0;
    for (b_i = 0; b_i <= k; b_i++) {
      if (cxy[b_i] < 0.0F) {
        trueCount++;
      }
    }
    r.set_size(trueCount);
    partialTrueCount = 0;
    for (b_i = 0; b_i <= k; b_i++) {
      if (cxy[b_i] < 0.0F) {
        r[partialTrueCount] = b_i + 1;
        partialTrueCount++;
      }
    }
    loop_ub = r.size(0);
    for (b_i = 0; b_i < loop_ub; b_i++) {
      cxy[r[b_i] - 1] = 0.0F;
    }
    siz_idx_0 = sigma * sigma;
    a = 1.5 * sigma;
    nx = I_45_45.size(0) * I_45_45.size(1);
    c45.set_size(I_45_45.size(0), I_45_45.size(1));
    if (static_cast<int>(nx < 3200)) {
      for (int j{0}; j < nx; j++) {
        c45[j] = std::abs(I_45_45[j]);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int j = 0; j < nx; j++) {
        c45[j] = std::abs(I_45_45[j]);
      }
    }
    nx = Ix.size(0) * Ix.size(1);
    c_y.set_size(Ix.size(0), Ix.size(1));
    if (static_cast<int>(nx < 3200)) {
      for (int j{0}; j < nx; j++) {
        c_y[j] = std::abs(Ix[j]);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int j = 0; j < nx; j++) {
        c_y[j] = std::abs(Ix[j]);
      }
    }
    nx = Iy.size(0) * Iy.size(1);
    checkerEdge.set_size(Iy.size(0), Iy.size(1));
    if (static_cast<int>(nx < 3200)) {
      for (int j{0}; j < nx; j++) {
        checkerEdge[j] = std::abs(Iy[j]);
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int j = 0; j < nx; j++) {
        checkerEdge[j] = std::abs(Iy[j]);
      }
    }
    if (c_y.size(0) == 1) {
      ny = checkerEdge.size(0);
    } else {
      ny = c_y.size(0);
    }
    if (c_y.size(1) == 1) {
      b_i = checkerEdge.size(1);
    } else {
      b_i = c_y.size(1);
    }
    if ((c_y.size(0) == checkerEdge.size(0)) &&
        (c_y.size(1) == checkerEdge.size(1)) && (c45.size(0) == ny) &&
        (c45.size(1) == b_i)) {
      loop_ub = c45.size(0) * c45.size(1);
      if (static_cast<int>(loop_ub < 3200)) {
        for (int j{0}; j < loop_ub; j++) {
          c45[j] = static_cast<float>(siz_idx_0) * c45[j] -
                   static_cast<float>(a) * (c_y[j] + checkerEdge[j]);
        }
      } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

        for (int j = 0; j < loop_ub; j++) {
          c45[j] = static_cast<float>(siz_idx_0) * c45[j] -
                   static_cast<float>(a) * (c_y[j] + checkerEdge[j]);
        }
      }
    } else {
      binary_expand_op(c45, siz_idx_0, a, c_y, checkerEdge);
    }
    k = c45.size(0) * c45.size(1) - 1;
    trueCount = 0;
    for (b_i = 0; b_i <= k; b_i++) {
      if (c45[b_i] < 0.0F) {
        trueCount++;
      }
    }
    r1.set_size(trueCount);
    partialTrueCount = 0;
    for (b_i = 0; b_i <= k; b_i++) {
      if (c45[b_i] < 0.0F) {
        r1[partialTrueCount] = b_i + 1;
        partialTrueCount++;
      }
    }
    loop_ub = r1.size(0);
    for (b_i = 0; b_i < loop_ub; b_i++) {
      c45[r1[b_i] - 1] = 0.0F;
    }
  } else {
    cell_wrap_11 derivFilters[12];
    int b_loop_ub;
    int b_nx;
    int c_loop_ub;
    int d_loop_ub;
    getDerivFilters(derivFilters);
    checkerboardPeaks.set_size(b_I.size(0), b_I.size(1));
    loop_ub = b_I.size(0) * b_I.size(1);
    if (static_cast<int>(loop_ub < 3200)) {
      for (int j{0}; j < loop_ub; j++) {
        checkerboardPeaks[j] = 0.0F;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int j = 0; j < loop_ub; j++) {
        checkerboardPeaks[j] = 0.0F;
      }
    }
    IfilterDot.set_size(Ig.size(0), Ig.size(1), 4);
    loop_ub = (Ig.size(0) * Ig.size(1)) << 2;
    if (static_cast<int>(loop_ub < 3200)) {
      for (int j{0}; j < loop_ub; j++) {
        IfilterDot[j] = 0.0F;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int j = 0; j < loop_ub; j++) {
        IfilterDot[j] = 0.0F;
      }
    }
    IfilterDdot.set_size(Ig.size(0), Ig.size(1), 4);
    loop_ub = (Ig.size(0) * Ig.size(1)) << 2;
    if (static_cast<int>(loop_ub < 3200)) {
      for (int j{0}; j < loop_ub; j++) {
        IfilterDdot[j] = 0.0F;
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int j = 0; j < loop_ub; j++) {
        IfilterDdot[j] = 0.0F;
      }
    }
    loop_ub = Ig.size(0) * Ig.size(1);
    b_loop_ub = Ig.size(0) * Ig.size(1);
    c_loop_ub = Ig.size(0) * Ig.size(1);
    d_loop_ub = Ig.size(0) * Ig.size(1);
    nx = Ix.size(0) * Ix.size(1);
    b_nx = Iy.size(0) * Iy.size(1);
    for (int filtIdx{0}; filtIdx < 3; filtIdx++) {
      float varargin_1;
      checkerEdge.set_size(Ig.size(0), Ig.size(1));
      for (b_i = 0; b_i < loop_ub; b_i++) {
        checkerEdge[b_i] = Ig[b_i];
      }
      imfilter(checkerEdge, derivFilters[filtIdx].f1);
      ny = checkerEdge.size(1);
      for (b_i = 0; b_i < ny; b_i++) {
        k = checkerEdge.size(0);
        for (partialTrueCount = 0; partialTrueCount < k; partialTrueCount++) {
          IfilterDot[partialTrueCount + IfilterDot.size(0) * b_i] =
              checkerEdge[partialTrueCount + checkerEdge.size(0) * b_i];
        }
      }
      checkerEdge.set_size(Ig.size(0), Ig.size(1));
      for (b_i = 0; b_i < b_loop_ub; b_i++) {
        checkerEdge[b_i] = Ig[b_i];
      }
      imfilter(checkerEdge, derivFilters[filtIdx + 3].f1);
      ny = checkerEdge.size(1);
      for (b_i = 0; b_i < ny; b_i++) {
        k = checkerEdge.size(0);
        for (partialTrueCount = 0; partialTrueCount < k; partialTrueCount++) {
          IfilterDot[(partialTrueCount + IfilterDot.size(0) * b_i) +
                     IfilterDot.size(0) * IfilterDot.size(1)] =
              checkerEdge[partialTrueCount + checkerEdge.size(0) * b_i];
        }
      }
      checkerEdge.set_size(Ig.size(0), Ig.size(1));
      for (b_i = 0; b_i < c_loop_ub; b_i++) {
        checkerEdge[b_i] = Ig[b_i];
      }
      imfilter(checkerEdge, derivFilters[filtIdx + 6].f1);
      ny = checkerEdge.size(1);
      for (b_i = 0; b_i < ny; b_i++) {
        k = checkerEdge.size(0);
        for (partialTrueCount = 0; partialTrueCount < k; partialTrueCount++) {
          IfilterDot[(partialTrueCount + IfilterDot.size(0) * b_i) +
                     IfilterDot.size(0) * IfilterDot.size(1) * 2] =
              checkerEdge[partialTrueCount + checkerEdge.size(0) * b_i];
        }
      }
      checkerEdge.set_size(Ig.size(0), Ig.size(1));
      for (b_i = 0; b_i < d_loop_ub; b_i++) {
        checkerEdge[b_i] = Ig[b_i];
      }
      imfilter(checkerEdge, derivFilters[filtIdx + 9].f1);
      ny = checkerEdge.size(1);
      for (b_i = 0; b_i < ny; b_i++) {
        k = checkerEdge.size(0);
        for (partialTrueCount = 0; partialTrueCount < k; partialTrueCount++) {
          IfilterDot[(partialTrueCount + IfilterDot.size(0) * b_i) +
                     IfilterDot.size(0) * IfilterDot.size(1) * 3] =
              checkerEdge[partialTrueCount + checkerEdge.size(0) * b_i];
        }
      }
      ny = IfilterDot.size(1);
      checkerEdge.set_size(IfilterDot.size(0), IfilterDot.size(1));
      for (b_i = 0; b_i < ny; b_i++) {
        k = IfilterDot.size(0);
        for (partialTrueCount = 0; partialTrueCount < k; partialTrueCount++) {
          checkerEdge[partialTrueCount + checkerEdge.size(0) * b_i] =
              IfilterDot[partialTrueCount + IfilterDot.size(0) * b_i];
        }
      }
      imfilter(checkerEdge, derivFilters[filtIdx].f1);
      ny = checkerEdge.size(1);
      for (b_i = 0; b_i < ny; b_i++) {
        k = checkerEdge.size(0);
        for (partialTrueCount = 0; partialTrueCount < k; partialTrueCount++) {
          IfilterDdot[partialTrueCount + IfilterDdot.size(0) * b_i] =
              checkerEdge[partialTrueCount + checkerEdge.size(0) * b_i];
        }
      }
      ny = IfilterDot.size(1);
      checkerEdge.set_size(IfilterDot.size(0), IfilterDot.size(1));
      for (b_i = 0; b_i < ny; b_i++) {
        k = IfilterDot.size(0);
        for (partialTrueCount = 0; partialTrueCount < k; partialTrueCount++) {
          checkerEdge[partialTrueCount + checkerEdge.size(0) * b_i] =
              IfilterDot[(partialTrueCount + IfilterDot.size(0) * b_i) +
                         IfilterDot.size(0) * IfilterDot.size(1)];
        }
      }
      imfilter(checkerEdge, derivFilters[filtIdx + 3].f1);
      ny = checkerEdge.size(1);
      for (b_i = 0; b_i < ny; b_i++) {
        k = checkerEdge.size(0);
        for (partialTrueCount = 0; partialTrueCount < k; partialTrueCount++) {
          IfilterDdot[(partialTrueCount + IfilterDdot.size(0) * b_i) +
                      IfilterDdot.size(0) * IfilterDdot.size(1)] =
              checkerEdge[partialTrueCount + checkerEdge.size(0) * b_i];
        }
      }
      ny = IfilterDot.size(1);
      checkerEdge.set_size(IfilterDot.size(0), IfilterDot.size(1));
      for (b_i = 0; b_i < ny; b_i++) {
        k = IfilterDot.size(0);
        for (partialTrueCount = 0; partialTrueCount < k; partialTrueCount++) {
          checkerEdge[partialTrueCount + checkerEdge.size(0) * b_i] =
              IfilterDot[partialTrueCount + IfilterDot.size(0) * b_i];
        }
      }
      imfilter(checkerEdge, derivFilters[filtIdx + 6].f1);
      ny = checkerEdge.size(1);
      for (b_i = 0; b_i < ny; b_i++) {
        k = checkerEdge.size(0);
        for (partialTrueCount = 0; partialTrueCount < k; partialTrueCount++) {
          IfilterDdot[(partialTrueCount + IfilterDdot.size(0) * b_i) +
                      IfilterDdot.size(0) * IfilterDdot.size(1) * 2] =
              checkerEdge[partialTrueCount + checkerEdge.size(0) * b_i];
        }
      }
      ny = IfilterDot.size(1);
      checkerEdge.set_size(IfilterDot.size(0), IfilterDot.size(1));
      for (b_i = 0; b_i < ny; b_i++) {
        k = IfilterDot.size(0);
        for (partialTrueCount = 0; partialTrueCount < k; partialTrueCount++) {
          checkerEdge[partialTrueCount + checkerEdge.size(0) * b_i] =
              IfilterDot[(partialTrueCount + IfilterDot.size(0) * b_i) +
                         IfilterDot.size(0) * IfilterDot.size(1)];
        }
      }
      imfilter(checkerEdge, derivFilters[filtIdx + 9].f1);
      ny = checkerEdge.size(1);
      for (b_i = 0; b_i < ny; b_i++) {
        k = checkerEdge.size(0);
        for (partialTrueCount = 0; partialTrueCount < k; partialTrueCount++) {
          IfilterDdot[(partialTrueCount + IfilterDdot.size(0) * b_i) +
                      IfilterDdot.size(0) * IfilterDdot.size(1) * 3] =
              checkerEdge[partialTrueCount + checkerEdge.size(0) * b_i];
        }
      }
      checkerEdgeComp.set_size(Ix.size(0), Ix.size(1));
      for (k = 0; k < nx; k++) {
        checkerEdgeComp[k] = std::abs(Ix[k]);
      }
      c_y.set_size(Iy.size(0), Iy.size(1));
      for (k = 0; k < b_nx; k++) {
        c_y[k] = std::abs(Iy[k]);
      }
      if ((checkerEdgeComp.size(0) == c_y.size(0)) &&
          (checkerEdgeComp.size(1) == c_y.size(1))) {
        ny = checkerEdgeComp.size(0) * checkerEdgeComp.size(1);
        for (b_i = 0; b_i < ny; b_i++) {
          checkerEdgeComp[b_i] = checkerEdgeComp[b_i] * 0.3F + c_y[b_i] * 0.3F;
        }
      } else {
        c_binary_expand_op(checkerEdgeComp, c_y);
      }
      ny = IfilterDdot.size(1);
      if ((IfilterDdot.size(0) == checkerEdgeComp.size(0)) &&
          (IfilterDdot.size(1) == checkerEdgeComp.size(1))) {
        checkerEdge.set_size(IfilterDdot.size(0), IfilterDdot.size(1));
        for (b_i = 0; b_i < ny; b_i++) {
          k = IfilterDdot.size(0);
          for (partialTrueCount = 0; partialTrueCount < k; partialTrueCount++) {
            checkerEdge[partialTrueCount + checkerEdge.size(0) * b_i] =
                (IfilterDdot[partialTrueCount + IfilterDdot.size(0) * b_i] +
                 IfilterDdot[(partialTrueCount + IfilterDdot.size(0) * b_i) +
                             IfilterDdot.size(0) * IfilterDdot.size(1)]) -
                checkerEdgeComp[partialTrueCount +
                                checkerEdgeComp.size(0) * b_i];
          }
        }
      } else {
        binary_expand_op(checkerEdge, IfilterDdot, checkerEdgeComp);
      }
      ny = IfilterDdot.size(1);
      if ((IfilterDdot.size(0) == checkerEdgeComp.size(0)) &&
          (IfilterDdot.size(1) == checkerEdgeComp.size(1))) {
        checkerEdgeComp.set_size(IfilterDdot.size(0), IfilterDdot.size(1));
        for (b_i = 0; b_i < ny; b_i++) {
          k = IfilterDdot.size(0);
          for (partialTrueCount = 0; partialTrueCount < k; partialTrueCount++) {
            checkerEdgeComp[partialTrueCount + checkerEdgeComp.size(0) * b_i] =
                (IfilterDdot[(partialTrueCount + IfilterDdot.size(0) * b_i) +
                             IfilterDdot.size(0) * IfilterDdot.size(1) * 2] +
                 IfilterDdot[(partialTrueCount + IfilterDdot.size(0) * b_i) +
                             IfilterDdot.size(0) * IfilterDdot.size(1) * 3]) -
                checkerEdgeComp[partialTrueCount +
                                checkerEdgeComp.size(0) * b_i];
          }
        }
      } else {
        binary_expand_op(checkerEdgeComp, IfilterDdot);
      }
      k = checkerEdge.size(0) * checkerEdge.size(1) - 1;
      trueCount = 0;
      for (b_i = 0; b_i <= k; b_i++) {
        if (checkerEdge[b_i] < 0.0F) {
          trueCount++;
        }
      }
      r.set_size(trueCount);
      partialTrueCount = 0;
      for (b_i = 0; b_i <= k; b_i++) {
        if (checkerEdge[b_i] < 0.0F) {
          r[partialTrueCount] = b_i + 1;
          partialTrueCount++;
        }
      }
      ny = r.size(0);
      for (b_i = 0; b_i < ny; b_i++) {
        checkerEdge[r[b_i] - 1] = 0.0F;
      }
      k = checkerEdgeComp.size(0) * checkerEdgeComp.size(1) - 1;
      trueCount = 0;
      for (b_i = 0; b_i <= k; b_i++) {
        if (checkerEdgeComp[b_i] < 0.0F) {
          trueCount++;
        }
      }
      r1.set_size(trueCount);
      partialTrueCount = 0;
      for (b_i = 0; b_i <= k; b_i++) {
        if (checkerEdgeComp[b_i] < 0.0F) {
          r1[partialTrueCount] = b_i + 1;
          partialTrueCount++;
        }
      }
      ny = r1.size(0);
      for (b_i = 0; b_i < ny; b_i++) {
        checkerEdgeComp[r1[b_i] - 1] = 0.0F;
      }
      c_y.set_size(checkerEdge.size(0), checkerEdge.size(1));
      ny = checkerEdge.size(0) * checkerEdge.size(1) - 1;
      for (b_i = 0; b_i <= ny; b_i++) {
        c_y[b_i] = checkerEdge[b_i];
      }
      imdilate(c_y, checkerEdge);
      c_y.set_size(checkerEdgeComp.size(0), checkerEdgeComp.size(1));
      ny = checkerEdgeComp.size(0) * checkerEdgeComp.size(1) - 1;
      for (b_i = 0; b_i <= ny; b_i++) {
        c_y[b_i] = checkerEdgeComp[b_i];
      }
      imdilate(c_y, checkerEdgeComp);
      ny = IfilterDdot.size(1);
      b_IfilterDdot.set_size(IfilterDdot.size(0), IfilterDdot.size(1), 2);
      for (b_i = 0; b_i < 2; b_i++) {
        for (partialTrueCount = 0; partialTrueCount < ny; partialTrueCount++) {
          k = IfilterDdot.size(0);
          for (trueCount = 0; trueCount < k; trueCount++) {
            b_IfilterDdot[(trueCount +
                           b_IfilterDdot.size(0) * partialTrueCount) +
                          b_IfilterDdot.size(0) * b_IfilterDdot.size(1) * b_i] =
                IfilterDdot[(trueCount +
                             IfilterDdot.size(0) * partialTrueCount) +
                            IfilterDdot.size(0) * IfilterDdot.size(1) * b_i];
          }
        }
      }
      bsxfun(b_IfilterDdot, checkerEdge, x);
      ny = x.size(0) * x.size(1) * 2;
      x.set_size(x.size(0), x.size(1), 2);
      for (b_i = 0; b_i < ny; b_i++) {
        varargin_1 = x[b_i];
        x[b_i] = std::fmax(varargin_1, 0.0F);
      }
      if ((x.size(0) == 0) || (x.size(1) == 0)) {
        c_y.set_size(x.size(0), x.size(1));
        ny = x.size(0) * x.size(1);
        for (b_i = 0; b_i < ny; b_i++) {
          c_y[b_i] = 0.0F;
        }
      } else {
        ny = x.size(0) * x.size(1);
        c_y.set_size(x.size(0), x.size(1));
        for (partialTrueCount = 0; partialTrueCount < ny; partialTrueCount++) {
          c_y[partialTrueCount] = x[partialTrueCount];
        }
        for (partialTrueCount = 0; partialTrueCount < ny; partialTrueCount++) {
          c_y[partialTrueCount] =
              c_y[partialTrueCount] + x[ny + partialTrueCount];
        }
      }
      ny = IfilterDdot.size(1);
      b_IfilterDdot.set_size(IfilterDdot.size(0), IfilterDdot.size(1), 2);
      for (b_i = 0; b_i < 2; b_i++) {
        for (partialTrueCount = 0; partialTrueCount < ny; partialTrueCount++) {
          k = IfilterDdot.size(0);
          for (trueCount = 0; trueCount < k; trueCount++) {
            b_IfilterDdot[(trueCount +
                           b_IfilterDdot.size(0) * partialTrueCount) +
                          b_IfilterDdot.size(0) * b_IfilterDdot.size(1) * b_i] =
                IfilterDdot[(trueCount +
                             IfilterDdot.size(0) * partialTrueCount) +
                            IfilterDdot.size(0) * IfilterDdot.size(1) *
                                (b_i + 2)];
          }
        }
      }
      bsxfun(b_IfilterDdot, checkerEdgeComp, x);
      ny = x.size(0) * x.size(1) * 2;
      x.set_size(x.size(0), x.size(1), 2);
      for (b_i = 0; b_i < ny; b_i++) {
        varargin_1 = x[b_i];
        x[b_i] = std::fmax(varargin_1, 0.0F);
      }
      if ((x.size(0) == 0) || (x.size(1) == 0)) {
        checkerEdge.set_size(x.size(0), x.size(1));
        ny = x.size(0) * x.size(1);
        for (b_i = 0; b_i < ny; b_i++) {
          checkerEdge[b_i] = 0.0F;
        }
      } else {
        ny = x.size(0) * x.size(1);
        checkerEdge.set_size(x.size(0), x.size(1));
        for (partialTrueCount = 0; partialTrueCount < ny; partialTrueCount++) {
          checkerEdge[partialTrueCount] = x[partialTrueCount];
        }
        for (partialTrueCount = 0; partialTrueCount < ny; partialTrueCount++) {
          checkerEdge[partialTrueCount] =
              checkerEdge[partialTrueCount] + x[ny + partialTrueCount];
        }
      }
      if (checkerboardPeaks.size(0) == 1) {
        ny = c_y.size(0);
      } else {
        ny = checkerboardPeaks.size(0);
      }
      if (checkerboardPeaks.size(1) == 1) {
        b_i = c_y.size(1);
      } else {
        b_i = checkerboardPeaks.size(1);
      }
      if ((checkerboardPeaks.size(0) == c_y.size(0)) &&
          (checkerboardPeaks.size(1) == c_y.size(1)) &&
          (ny == checkerEdge.size(0)) && (b_i == checkerEdge.size(1))) {
        ny = checkerboardPeaks.size(0) * checkerboardPeaks.size(1);
        for (b_i = 0; b_i < ny; b_i++) {
          checkerboardPeaks[b_i] =
              (checkerboardPeaks[b_i] + c_y[b_i]) + checkerEdge[b_i];
        }
      } else {
        b_binary_expand_op(checkerboardPeaks, c_y, checkerEdge);
      }
    }
    imerode(checkerboardPeaks, cxy);
  }
}

} // namespace checkerboard
} // namespace calibration
} // namespace internal
} // namespace vision
} // namespace coder

//
// File trailer for secondDerivCornerMetric.cpp
//
// [EOF]
//
