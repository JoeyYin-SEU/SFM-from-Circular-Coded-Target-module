//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: algbwmorph.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "algbwmorph.h"
#include "get_chessborad_pixel_rtwutil.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "libmwbwlookup_tbb.h"
#include "omp.h"

// Function Declarations
static void g_binary_expand_op(coder::array<boolean_T, 2U> &in1,
                               const coder::array<boolean_T, 2U> &in2);

// Function Definitions
//
// Arguments    : coder::array<boolean_T, 2U> &in1
//                const coder::array<boolean_T, 2U> &in2
// Return Type  : void
//
static void g_binary_expand_op(coder::array<boolean_T, 2U> &in1,
                               const coder::array<boolean_T, 2U> &in2)
{
  coder::array<boolean_T, 2U> b_in2;
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
          (in2[i1 * stride_0_0 + in2.size(0) * aux_0_1] &&
           (!in1[i1 * stride_1_0 + in1.size(0) * aux_1_1]));
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
// Arguments    : ::coder::array<boolean_T, 2U> &bw
// Return Type  : void
//
namespace coder {
namespace images {
namespace internal {
void bwmorphApplyOnce(::coder::array<boolean_T, 2U> &bw)
{
  static const boolean_T lut[512]{
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, true,  true,  true,
      true,  false, true,  true,  true,  true,  true,  true,  false, false,
      true,  true,  false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, true,  false,
      true,  true,  true,  false, true,  true,  false, false, true,  true,
      false, false, true,  true,  false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      true,  false, false, false, false, false, false, false, true,  true,
      true,  true,  false, false, true,  true,  false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, true,  true,  false, false, true,  true,  false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, true,  false, false, false, false, false,
      false, false, true,  true,  true,  true,  false, false, true,  true,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, true,  false, true,  true,
      true,  false, true,  true,  true,  true,  false, false, true,  true,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, true,  false,
      false, false, false, false, false, false, true,  true,  true,  true,
      false, false, true,  true,  false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      true,  false, true,  true,  true,  false, true,  true,  true,  true,
      false, false, true,  true,  false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, true,  false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, true,  false, true,  true,  true,  false,
      true,  true,  false, false, true,  true,  false, false, true,  true,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, true,  true,
      false, false, true,  true,  false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      true,  false, false, false, false, false, false, false, true,  true,
      true,  true,  false, false, true,  true,  false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, true,  false, true,  true,  true,  false, true,  true,
      true,  true,  false, false, true,  true,  false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, true,  false, false, false, false, false,
      false, false, true,  true,  true,  true,  false, false, true,  true,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, true,  false, true,  true,
      true,  false, true,  true,  true,  true,  false, false, true,  true,
      false, false};
  array<boolean_T, 2U> m;
  double sizeIn[2];
  int b_loop_ub;
  int i;
  int i10;
  int i2;
  int i3;
  int i4;
  int i5;
  int i6;
  int i7;
  int i8;
  int i9;
  int loop_ub;
  m.set_size(bw.size(0), bw.size(1));
  if ((bw.size(0) != 0) && (bw.size(1) != 0)) {
    sizeIn[0] = bw.size(0);
    sizeIn[1] = bw.size(1);
    bwlookup_tbb_boolean(&bw[0], &sizeIn[0], 2.0, &lut[0], 512.0, &m[0]);
  }
  if ((bw.size(0) == m.size(0)) && (bw.size(1) == m.size(1))) {
    loop_ub = bw.size(0) * bw.size(1);
    m.set_size(bw.size(0), bw.size(1));
    if (static_cast<int>(loop_ub < 3200)) {
      for (int i1{0}; i1 < loop_ub; i1++) {
        m[i1] = (bw[i1] && (!m[i1]));
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i1 = 0; i1 < loop_ub; i1++) {
        m[i1] = (bw[i1] && (!m[i1]));
      }
    }
  } else {
    g_binary_expand_op(m, bw);
  }
  if (m.size(0) < 1) {
    i = 1;
    i2 = 0;
  } else {
    i = 2;
    i2 = m.size(0);
  }
  if (m.size(1) < 1) {
    i3 = 1;
    i4 = 0;
  } else {
    i3 = 2;
    i4 = m.size(1);
  }
  if (bw.size(0) < 1) {
    i5 = 1;
  } else {
    i5 = 2;
  }
  if (bw.size(1) < 1) {
    i6 = 1;
  } else {
    i6 = 2;
  }
  loop_ub = div_s32(i4 - 1, i3);
  for (i4 = 0; i4 <= loop_ub; i4++) {
    b_loop_ub = div_s32(i2 - 1, i);
    for (i7 = 0; i7 <= b_loop_ub; i7++) {
      bw[i5 * i7 + bw.size(0) * (i6 * i4)] = m[i * i7 + m.size(0) * (i3 * i4)];
    }
  }
  m.set_size(bw.size(0), bw.size(1));
  if ((bw.size(0) != 0) && (bw.size(1) != 0)) {
    sizeIn[0] = bw.size(0);
    sizeIn[1] = bw.size(1);
    bwlookup_tbb_boolean(&bw[0], &sizeIn[0], 2.0, &lut[0], 512.0, &m[0]);
  }
  if ((bw.size(0) == m.size(0)) && (bw.size(1) == m.size(1))) {
    loop_ub = bw.size(0) * bw.size(1);
    m.set_size(bw.size(0), bw.size(1));
    if (static_cast<int>(loop_ub < 3200)) {
      for (int i1{0}; i1 < loop_ub; i1++) {
        m[i1] = (bw[i1] && (!m[i1]));
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i1 = 0; i1 < loop_ub; i1++) {
        m[i1] = (bw[i1] && (!m[i1]));
      }
    }
  } else {
    g_binary_expand_op(m, bw);
  }
  if (m.size(0) < 2) {
    i = 0;
    i2 = 1;
    i3 = 0;
  } else {
    i = 1;
    i2 = 2;
    i3 = m.size(0);
  }
  if (m.size(1) < 2) {
    i4 = 0;
    i5 = 1;
    i6 = 0;
  } else {
    i4 = 1;
    i5 = 2;
    i6 = m.size(1);
  }
  if (bw.size(0) < 2) {
    i7 = 0;
    i8 = 1;
  } else {
    i7 = 1;
    i8 = 2;
  }
  if (bw.size(1) < 2) {
    i9 = 0;
    i10 = 1;
  } else {
    i9 = 1;
    i10 = 2;
  }
  loop_ub = div_s32((i6 - i4) - 1, i5);
  for (i6 = 0; i6 <= loop_ub; i6++) {
    b_loop_ub = div_s32((i3 - i) - 1, i2);
    for (int i11{0}; i11 <= b_loop_ub; i11++) {
      bw[(i7 + i8 * i11) + bw.size(0) * (i9 + i10 * i6)] =
          m[(i + i2 * i11) + m.size(0) * (i4 + i5 * i6)];
    }
  }
  m.set_size(bw.size(0), bw.size(1));
  if ((bw.size(0) != 0) && (bw.size(1) != 0)) {
    sizeIn[0] = bw.size(0);
    sizeIn[1] = bw.size(1);
    bwlookup_tbb_boolean(&bw[0], &sizeIn[0], 2.0, &lut[0], 512.0, &m[0]);
  }
  if ((bw.size(0) == m.size(0)) && (bw.size(1) == m.size(1))) {
    loop_ub = bw.size(0) * bw.size(1);
    m.set_size(bw.size(0), bw.size(1));
    if (static_cast<int>(loop_ub < 3200)) {
      for (int i1{0}; i1 < loop_ub; i1++) {
        m[i1] = (bw[i1] && (!m[i1]));
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i1 = 0; i1 < loop_ub; i1++) {
        m[i1] = (bw[i1] && (!m[i1]));
      }
    }
  } else {
    g_binary_expand_op(m, bw);
  }
  if (m.size(0) < 1) {
    i = 1;
    i2 = 0;
  } else {
    i = 2;
    i2 = m.size(0);
  }
  if (m.size(1) < 2) {
    i3 = 0;
    i4 = 1;
    i5 = 0;
  } else {
    i3 = 1;
    i4 = 2;
    i5 = m.size(1);
  }
  if (bw.size(0) < 1) {
    i6 = 1;
  } else {
    i6 = 2;
  }
  if (bw.size(1) < 2) {
    i7 = 0;
    i8 = 1;
  } else {
    i7 = 1;
    i8 = 2;
  }
  loop_ub = div_s32((i5 - i3) - 1, i4);
  for (i5 = 0; i5 <= loop_ub; i5++) {
    b_loop_ub = div_s32(i2 - 1, i);
    for (i9 = 0; i9 <= b_loop_ub; i9++) {
      bw[i6 * i9 + bw.size(0) * (i7 + i8 * i5)] =
          m[i * i9 + m.size(0) * (i3 + i4 * i5)];
    }
  }
  m.set_size(bw.size(0), bw.size(1));
  if ((bw.size(0) != 0) && (bw.size(1) != 0)) {
    sizeIn[0] = bw.size(0);
    sizeIn[1] = bw.size(1);
    bwlookup_tbb_boolean(&bw[0], &sizeIn[0], 2.0, &lut[0], 512.0, &m[0]);
  }
  if ((bw.size(0) == m.size(0)) && (bw.size(1) == m.size(1))) {
    loop_ub = bw.size(0) * bw.size(1);
    m.set_size(bw.size(0), bw.size(1));
    if (static_cast<int>(loop_ub < 3200)) {
      for (int i1{0}; i1 < loop_ub; i1++) {
        m[i1] = (bw[i1] && (!m[i1]));
      }
    } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

      for (int i1 = 0; i1 < loop_ub; i1++) {
        m[i1] = (bw[i1] && (!m[i1]));
      }
    }
  } else {
    g_binary_expand_op(m, bw);
  }
  if (m.size(0) < 2) {
    i = 0;
    i2 = 1;
    i3 = 0;
  } else {
    i = 1;
    i2 = 2;
    i3 = m.size(0);
  }
  if (m.size(1) < 1) {
    i4 = 1;
    i5 = 0;
  } else {
    i4 = 2;
    i5 = m.size(1);
  }
  if (bw.size(0) < 2) {
    i6 = 0;
    i7 = 1;
  } else {
    i6 = 1;
    i7 = 2;
  }
  if (bw.size(1) < 1) {
    i8 = 1;
  } else {
    i8 = 2;
  }
  loop_ub = div_s32(i5 - 1, i4);
  for (i5 = 0; i5 <= loop_ub; i5++) {
    b_loop_ub = div_s32((i3 - i) - 1, i2);
    for (i9 = 0; i9 <= b_loop_ub; i9++) {
      bw[(i6 + i7 * i9) + bw.size(0) * (i8 * i5)] =
          m[(i + i2 * i9) + m.size(0) * (i4 * i5)];
    }
  }
}

} // namespace internal
} // namespace images
} // namespace coder

//
// File trailer for algbwmorph.cpp
//
// [EOF]
//
