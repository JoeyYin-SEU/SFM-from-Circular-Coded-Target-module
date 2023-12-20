//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: sort.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "sort.h"
#include "rt_nonfinite.h"
#include "sortIdx.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
//
// Arguments    : ::coder::array<float, 1U> &x
//                ::coder::array<int, 1U> &idx
// Return Type  : void
//
namespace coder {
namespace internal {
void sort(::coder::array<float, 1U> &x, ::coder::array<int, 1U> &idx)
{
  array<float, 1U> vwork;
  array<float, 1U> xwork;
  array<int, 1U> iidx;
  array<int, 1U> iwork;
  int dim;
  int i;
  int vlen;
  int vstride;
  dim = 0;
  if (x.size(0) != 1) {
    dim = -1;
  }
  if (dim + 2 <= 1) {
    i = x.size(0);
  } else {
    i = 1;
  }
  vlen = i - 1;
  vwork.set_size(i);
  idx.set_size(x.size(0));
  vstride = 1;
  for (int k{0}; k <= dim; k++) {
    vstride *= x.size(0);
  }
  for (int j{0}; j < vstride; j++) {
    for (int k{0}; k <= vlen; k++) {
      vwork[k] = x[j + k * vstride];
    }
    iidx.set_size(vwork.size(0));
    dim = vwork.size(0);
    for (i = 0; i < dim; i++) {
      iidx[i] = 0;
    }
    if (vwork.size(0) != 0) {
      float x4[4];
      int idx4[4];
      int i1;
      int i2;
      int i3;
      int i4;
      int iidx_tmp;
      int n;
      int nNaNs;
      int nNonNaN;
      n = vwork.size(0);
      x4[0] = 0.0F;
      idx4[0] = 0;
      x4[1] = 0.0F;
      idx4[1] = 0;
      x4[2] = 0.0F;
      idx4[2] = 0;
      x4[3] = 0.0F;
      idx4[3] = 0;
      iwork.set_size(vwork.size(0));
      dim = vwork.size(0);
      for (i = 0; i < dim; i++) {
        iwork[i] = 0;
      }
      xwork.set_size(vwork.size(0));
      dim = vwork.size(0);
      for (i = 0; i < dim; i++) {
        xwork[i] = 0.0F;
      }
      nNaNs = 0;
      dim = 0;
      for (int k{0}; k < n; k++) {
        if (std::isnan(vwork[k])) {
          iidx_tmp = (n - nNaNs) - 1;
          iidx[iidx_tmp] = k + 1;
          xwork[iidx_tmp] = vwork[k];
          nNaNs++;
        } else {
          dim++;
          idx4[dim - 1] = k + 1;
          x4[dim - 1] = vwork[k];
          if (dim == 4) {
            float f;
            float f1;
            signed char b_i1;
            signed char b_i2;
            signed char b_i3;
            signed char b_i4;
            dim = k - nNaNs;
            if (x4[0] >= x4[1]) {
              i1 = 1;
              i2 = 2;
            } else {
              i1 = 2;
              i2 = 1;
            }
            if (x4[2] >= x4[3]) {
              i3 = 3;
              i4 = 4;
            } else {
              i3 = 4;
              i4 = 3;
            }
            f = x4[i1 - 1];
            f1 = x4[i3 - 1];
            if (f >= f1) {
              f = x4[i2 - 1];
              if (f >= f1) {
                b_i1 = static_cast<signed char>(i1);
                b_i2 = static_cast<signed char>(i2);
                b_i3 = static_cast<signed char>(i3);
                b_i4 = static_cast<signed char>(i4);
              } else if (f >= x4[i4 - 1]) {
                b_i1 = static_cast<signed char>(i1);
                b_i2 = static_cast<signed char>(i3);
                b_i3 = static_cast<signed char>(i2);
                b_i4 = static_cast<signed char>(i4);
              } else {
                b_i1 = static_cast<signed char>(i1);
                b_i2 = static_cast<signed char>(i3);
                b_i3 = static_cast<signed char>(i4);
                b_i4 = static_cast<signed char>(i2);
              }
            } else {
              f1 = x4[i4 - 1];
              if (f >= f1) {
                if (x4[i2 - 1] >= f1) {
                  b_i1 = static_cast<signed char>(i3);
                  b_i2 = static_cast<signed char>(i1);
                  b_i3 = static_cast<signed char>(i2);
                  b_i4 = static_cast<signed char>(i4);
                } else {
                  b_i1 = static_cast<signed char>(i3);
                  b_i2 = static_cast<signed char>(i1);
                  b_i3 = static_cast<signed char>(i4);
                  b_i4 = static_cast<signed char>(i2);
                }
              } else {
                b_i1 = static_cast<signed char>(i3);
                b_i2 = static_cast<signed char>(i4);
                b_i3 = static_cast<signed char>(i1);
                b_i4 = static_cast<signed char>(i2);
              }
            }
            iidx[dim - 3] = idx4[b_i1 - 1];
            iidx[dim - 2] = idx4[b_i2 - 1];
            iidx[dim - 1] = idx4[b_i3 - 1];
            iidx[dim] = idx4[b_i4 - 1];
            vwork[dim - 3] = x4[b_i1 - 1];
            vwork[dim - 2] = x4[b_i2 - 1];
            vwork[dim - 1] = x4[b_i3 - 1];
            vwork[dim] = x4[b_i4 - 1];
            dim = 0;
          }
        }
      }
      i3 = vwork.size(0) - nNaNs;
      if (dim > 0) {
        signed char perm[4];
        perm[1] = 0;
        perm[2] = 0;
        perm[3] = 0;
        if (dim == 1) {
          perm[0] = 1;
        } else if (dim == 2) {
          if (x4[0] >= x4[1]) {
            perm[0] = 1;
            perm[1] = 2;
          } else {
            perm[0] = 2;
            perm[1] = 1;
          }
        } else if (x4[0] >= x4[1]) {
          if (x4[1] >= x4[2]) {
            perm[0] = 1;
            perm[1] = 2;
            perm[2] = 3;
          } else if (x4[0] >= x4[2]) {
            perm[0] = 1;
            perm[1] = 3;
            perm[2] = 2;
          } else {
            perm[0] = 3;
            perm[1] = 1;
            perm[2] = 2;
          }
        } else if (x4[0] >= x4[2]) {
          perm[0] = 2;
          perm[1] = 1;
          perm[2] = 3;
        } else if (x4[1] >= x4[2]) {
          perm[0] = 2;
          perm[1] = 3;
          perm[2] = 1;
        } else {
          perm[0] = 3;
          perm[1] = 2;
          perm[2] = 1;
        }
        i = static_cast<unsigned char>(dim);
        for (int k{0}; k < i; k++) {
          iidx_tmp = perm[k] - 1;
          i1 = (i3 - dim) + k;
          iidx[i1] = idx4[iidx_tmp];
          vwork[i1] = x4[iidx_tmp];
        }
      }
      dim = nNaNs >> 1;
      for (int k{0}; k < dim; k++) {
        i1 = i3 + k;
        i2 = iidx[i1];
        iidx_tmp = (n - k) - 1;
        iidx[i1] = iidx[iidx_tmp];
        iidx[iidx_tmp] = i2;
        vwork[i1] = xwork[iidx_tmp];
        vwork[iidx_tmp] = xwork[i1];
      }
      if ((nNaNs & 1) != 0) {
        dim += i3;
        vwork[dim] = xwork[dim];
      }
      nNonNaN = vwork.size(0) - nNaNs;
      dim = 2;
      if (nNonNaN > 1) {
        if (vwork.size(0) >= 256) {
          int nBlocks;
          nBlocks = nNonNaN >> 8;
          if (nBlocks > 0) {
            for (int b{0}; b < nBlocks; b++) {
              float b_xwork[256];
              int b_iwork[256];
              i4 = (b << 8) - 1;
              for (int b_b{0}; b_b < 6; b_b++) {
                int bLen2;
                n = 1 << (b_b + 2);
                bLen2 = n << 1;
                i = 256 >> (b_b + 3);
                for (int k{0}; k < i; k++) {
                  i2 = (i4 + k * bLen2) + 1;
                  for (i1 = 0; i1 < bLen2; i1++) {
                    dim = i2 + i1;
                    b_iwork[i1] = iidx[dim];
                    b_xwork[i1] = vwork[dim];
                  }
                  i3 = 0;
                  i1 = n;
                  dim = i2 - 1;
                  int exitg1;
                  do {
                    exitg1 = 0;
                    dim++;
                    if (b_xwork[i3] >= b_xwork[i1]) {
                      iidx[dim] = b_iwork[i3];
                      vwork[dim] = b_xwork[i3];
                      if (i3 + 1 < n) {
                        i3++;
                      } else {
                        exitg1 = 1;
                      }
                    } else {
                      iidx[dim] = b_iwork[i1];
                      vwork[dim] = b_xwork[i1];
                      if (i1 + 1 < bLen2) {
                        i1++;
                      } else {
                        dim -= i3;
                        for (i1 = i3 + 1; i1 <= n; i1++) {
                          iidx_tmp = dim + i1;
                          iidx[iidx_tmp] = b_iwork[i1 - 1];
                          vwork[iidx_tmp] = b_xwork[i1 - 1];
                        }
                        exitg1 = 1;
                      }
                    }
                  } while (exitg1 == 0);
                }
              }
            }
            dim = nBlocks << 8;
            i1 = nNonNaN - dim;
            if (i1 > 0) {
              merge_block(iidx, vwork, dim, i1, 2, iwork, xwork);
            }
            dim = 8;
          }
        }
        merge_block(iidx, vwork, 0, nNonNaN, dim, iwork, xwork);
      }
      if ((nNaNs > 0) && (nNonNaN > 0)) {
        for (int k{0}; k < nNaNs; k++) {
          dim = nNonNaN + k;
          xwork[k] = vwork[dim];
          iwork[k] = iidx[dim];
        }
        for (int k{nNonNaN}; k >= 1; k--) {
          dim = (nNaNs + k) - 1;
          vwork[dim] = vwork[k - 1];
          iidx[dim] = iidx[k - 1];
        }
        for (int k{0}; k < nNaNs; k++) {
          vwork[k] = xwork[k];
          iidx[k] = iwork[k];
        }
      }
    }
    for (int k{0}; k <= vlen; k++) {
      i = j + k * vstride;
      x[i] = vwork[k];
      idx[i] = iidx[k];
    }
  }
}

//
// Arguments    : ::coder::array<int, 1U> &x
//                ::coder::array<int, 1U> &idx
// Return Type  : void
//
void sort(::coder::array<int, 1U> &x, ::coder::array<int, 1U> &idx)
{
  array<int, 1U> b_iwork;
  array<int, 1U> iidx;
  array<int, 1U> iwork;
  array<int, 1U> vwork;
  array<int, 1U> xwork;
  int dim;
  int i;
  int vlen;
  int vstride;
  dim = 0;
  if (x.size(0) != 1) {
    dim = -1;
  }
  if (dim + 2 <= 1) {
    i = x.size(0);
  } else {
    i = 1;
  }
  vlen = i - 1;
  vwork.set_size(i);
  idx.set_size(x.size(0));
  vstride = 1;
  for (int k{0}; k <= dim; k++) {
    vstride *= x.size(0);
  }
  for (int j{0}; j < vstride; j++) {
    for (int k{0}; k <= vlen; k++) {
      vwork[k] = x[j + k * vstride];
    }
    iidx.set_size(vwork.size(0));
    dim = vwork.size(0);
    for (i = 0; i < dim; i++) {
      iidx[i] = 0;
    }
    if (vwork.size(0) != 0) {
      int idx4[4];
      int x4[4];
      int b_i;
      int i1;
      int i2;
      int i4;
      int nLeft;
      int nQuartets;
      x4[0] = 0;
      idx4[0] = 0;
      x4[1] = 0;
      idx4[1] = 0;
      x4[2] = 0;
      idx4[2] = 0;
      x4[3] = 0;
      idx4[3] = 0;
      iwork.set_size(vwork.size(0));
      dim = vwork.size(0);
      for (i = 0; i < dim; i++) {
        iwork[i] = 0;
      }
      xwork.set_size(vwork.size(0));
      dim = vwork.size(0);
      for (i = 0; i < dim; i++) {
        xwork[i] = 0;
      }
      nQuartets = vwork.size(0) >> 2;
      for (int b_j{0}; b_j < nQuartets; b_j++) {
        signed char b_i1;
        signed char b_i2;
        signed char b_i4;
        signed char i3;
        b_i = b_j << 2;
        idx4[0] = b_i + 1;
        idx4[1] = b_i + 2;
        idx4[2] = b_i + 3;
        idx4[3] = b_i + 4;
        x4[0] = vwork[b_i];
        dim = vwork[b_i + 1];
        x4[1] = dim;
        i4 = vwork[b_i + 2];
        x4[2] = i4;
        nLeft = vwork[b_i + 3];
        x4[3] = nLeft;
        if (vwork[b_i] <= dim) {
          i1 = 1;
          i2 = 2;
        } else {
          i1 = 2;
          i2 = 1;
        }
        if (i4 <= nLeft) {
          dim = 3;
          i4 = 4;
        } else {
          dim = 4;
          i4 = 3;
        }
        i = x4[i1 - 1];
        nLeft = x4[dim - 1];
        if (i <= nLeft) {
          i = x4[i2 - 1];
          if (i <= nLeft) {
            b_i1 = static_cast<signed char>(i1);
            b_i2 = static_cast<signed char>(i2);
            i3 = static_cast<signed char>(dim);
            b_i4 = static_cast<signed char>(i4);
          } else if (i <= x4[i4 - 1]) {
            b_i1 = static_cast<signed char>(i1);
            b_i2 = static_cast<signed char>(dim);
            i3 = static_cast<signed char>(i2);
            b_i4 = static_cast<signed char>(i4);
          } else {
            b_i1 = static_cast<signed char>(i1);
            b_i2 = static_cast<signed char>(dim);
            i3 = static_cast<signed char>(i4);
            b_i4 = static_cast<signed char>(i2);
          }
        } else {
          nLeft = x4[i4 - 1];
          if (i <= nLeft) {
            if (x4[i2 - 1] <= nLeft) {
              b_i1 = static_cast<signed char>(dim);
              b_i2 = static_cast<signed char>(i1);
              i3 = static_cast<signed char>(i2);
              b_i4 = static_cast<signed char>(i4);
            } else {
              b_i1 = static_cast<signed char>(dim);
              b_i2 = static_cast<signed char>(i1);
              i3 = static_cast<signed char>(i4);
              b_i4 = static_cast<signed char>(i2);
            }
          } else {
            b_i1 = static_cast<signed char>(dim);
            b_i2 = static_cast<signed char>(i4);
            i3 = static_cast<signed char>(i1);
            b_i4 = static_cast<signed char>(i2);
          }
        }
        iidx[b_i] = idx4[b_i1 - 1];
        iidx[b_i + 1] = idx4[b_i2 - 1];
        iidx[b_i + 2] = idx4[i3 - 1];
        iidx[b_i + 3] = idx4[b_i4 - 1];
        vwork[b_i] = x4[b_i1 - 1];
        vwork[b_i + 1] = x4[b_i2 - 1];
        vwork[b_i + 2] = x4[i3 - 1];
        vwork[b_i + 3] = x4[b_i4 - 1];
      }
      i4 = nQuartets << 2;
      nLeft = vwork.size(0) - i4;
      if (nLeft > 0) {
        signed char perm[4];
        for (int k{0}; k < nLeft; k++) {
          dim = i4 + k;
          idx4[k] = dim + 1;
          x4[k] = vwork[dim];
        }
        perm[1] = 0;
        perm[2] = 0;
        perm[3] = 0;
        if (nLeft == 1) {
          perm[0] = 1;
        } else if (nLeft == 2) {
          if (x4[0] <= x4[1]) {
            perm[0] = 1;
            perm[1] = 2;
          } else {
            perm[0] = 2;
            perm[1] = 1;
          }
        } else if (x4[0] <= x4[1]) {
          if (x4[1] <= x4[2]) {
            perm[0] = 1;
            perm[1] = 2;
            perm[2] = 3;
          } else if (x4[0] <= x4[2]) {
            perm[0] = 1;
            perm[1] = 3;
            perm[2] = 2;
          } else {
            perm[0] = 3;
            perm[1] = 1;
            perm[2] = 2;
          }
        } else if (x4[0] <= x4[2]) {
          perm[0] = 2;
          perm[1] = 1;
          perm[2] = 3;
        } else if (x4[1] <= x4[2]) {
          perm[0] = 2;
          perm[1] = 3;
          perm[2] = 1;
        } else {
          perm[0] = 3;
          perm[1] = 2;
          perm[2] = 1;
        }
        for (int k{0}; k < nLeft; k++) {
          i1 = perm[k] - 1;
          dim = i4 + k;
          iidx[dim] = idx4[i1];
          vwork[dim] = x4[i1];
        }
      }
      i4 = 2;
      if (vwork.size(0) > 1) {
        if (vwork.size(0) >= 256) {
          nQuartets = vwork.size(0) >> 8;
          for (int b{0}; b < nQuartets; b++) {
            int b_xwork[256];
            int c_iwork[256];
            b_i = (b << 8) - 1;
            for (int b_b{0}; b_b < 6; b_b++) {
              int bLen;
              int bLen2;
              bLen = 1 << (b_b + 2);
              bLen2 = bLen << 1;
              i = 256 >> (b_b + 3);
              for (int k{0}; k < i; k++) {
                i4 = (b_i + k * bLen2) + 1;
                for (int b_j{0}; b_j < bLen2; b_j++) {
                  dim = i4 + b_j;
                  c_iwork[b_j] = iidx[dim];
                  b_xwork[b_j] = vwork[dim];
                }
                i2 = 0;
                nLeft = bLen;
                dim = i4 - 1;
                int exitg1;
                do {
                  exitg1 = 0;
                  dim++;
                  if (b_xwork[i2] <= b_xwork[nLeft]) {
                    iidx[dim] = c_iwork[i2];
                    vwork[dim] = b_xwork[i2];
                    if (i2 + 1 < bLen) {
                      i2++;
                    } else {
                      exitg1 = 1;
                    }
                  } else {
                    iidx[dim] = c_iwork[nLeft];
                    vwork[dim] = b_xwork[nLeft];
                    if (nLeft + 1 < bLen2) {
                      nLeft++;
                    } else {
                      dim -= i2;
                      for (int b_j{i2 + 1}; b_j <= bLen; b_j++) {
                        i1 = dim + b_j;
                        iidx[i1] = c_iwork[b_j - 1];
                        vwork[i1] = b_xwork[b_j - 1];
                      }
                      exitg1 = 1;
                    }
                  }
                } while (exitg1 == 0);
              }
            }
          }
          dim = nQuartets << 8;
          i4 = vwork.size(0) - dim;
          if (i4 > 0) {
            merge_block(iidx, vwork, dim, i4, 2, iwork, xwork);
          }
          i4 = 8;
        }
        dim = iwork.size(0);
        b_iwork.set_size(iwork.size(0));
        for (i = 0; i < dim; i++) {
          b_iwork[i] = iwork[i];
        }
        dim = xwork.size(0);
        iwork.set_size(xwork.size(0));
        for (i = 0; i < dim; i++) {
          iwork[i] = xwork[i];
        }
        merge_block(iidx, vwork, 0, vwork.size(0), i4, b_iwork, iwork);
      }
    }
    for (int k{0}; k <= vlen; k++) {
      i = j + k * vstride;
      x[i] = vwork[k];
      idx[i] = iidx[k];
    }
  }
}

//
// Arguments    : ::coder::array<double, 1U> &x
//                ::coder::array<int, 1U> &idx
// Return Type  : void
//
void sort(::coder::array<double, 1U> &x, ::coder::array<int, 1U> &idx)
{
  array<double, 1U> b_xwork;
  array<double, 1U> vwork;
  array<double, 1U> xwork;
  array<int, 1U> b_iwork;
  array<int, 1U> iidx;
  array<int, 1U> iwork;
  int dim;
  int i;
  int vlen;
  int vstride;
  dim = 0;
  if (x.size(0) != 1) {
    dim = -1;
  }
  if (dim + 2 <= 1) {
    i = x.size(0);
  } else {
    i = 1;
  }
  vlen = i - 1;
  vwork.set_size(i);
  idx.set_size(x.size(0));
  vstride = 1;
  for (int k{0}; k <= dim; k++) {
    vstride *= x.size(0);
  }
  for (int j{0}; j < vstride; j++) {
    for (int k{0}; k <= vlen; k++) {
      vwork[k] = x[j + k * vstride];
    }
    iidx.set_size(vwork.size(0));
    dim = vwork.size(0);
    for (i = 0; i < dim; i++) {
      iidx[i] = 0;
    }
    if (vwork.size(0) != 0) {
      double x4[4];
      int idx4[4];
      int bLen2;
      int i1;
      int i2;
      int i3;
      int i4;
      int iidx_tmp;
      int n;
      int nNonNaN;
      n = vwork.size(0);
      x4[0] = 0.0;
      idx4[0] = 0;
      x4[1] = 0.0;
      idx4[1] = 0;
      x4[2] = 0.0;
      idx4[2] = 0;
      x4[3] = 0.0;
      idx4[3] = 0;
      iwork.set_size(vwork.size(0));
      dim = vwork.size(0);
      for (i = 0; i < dim; i++) {
        iwork[i] = 0;
      }
      xwork.set_size(vwork.size(0));
      dim = vwork.size(0);
      for (i = 0; i < dim; i++) {
        xwork[i] = 0.0;
      }
      bLen2 = 0;
      dim = 0;
      for (int k{0}; k < n; k++) {
        if (std::isnan(vwork[k])) {
          iidx_tmp = (n - bLen2) - 1;
          iidx[iidx_tmp] = k + 1;
          xwork[iidx_tmp] = vwork[k];
          bLen2++;
        } else {
          dim++;
          idx4[dim - 1] = k + 1;
          x4[dim - 1] = vwork[k];
          if (dim == 4) {
            double d;
            double d1;
            signed char b_i1;
            signed char b_i2;
            signed char b_i3;
            signed char b_i4;
            dim = k - bLen2;
            if (x4[0] <= x4[1]) {
              i1 = 1;
              i2 = 2;
            } else {
              i1 = 2;
              i2 = 1;
            }
            if (x4[2] <= x4[3]) {
              i3 = 3;
              i4 = 4;
            } else {
              i3 = 4;
              i4 = 3;
            }
            d = x4[i1 - 1];
            d1 = x4[i3 - 1];
            if (d <= d1) {
              d = x4[i2 - 1];
              if (d <= d1) {
                b_i1 = static_cast<signed char>(i1);
                b_i2 = static_cast<signed char>(i2);
                b_i3 = static_cast<signed char>(i3);
                b_i4 = static_cast<signed char>(i4);
              } else if (d <= x4[i4 - 1]) {
                b_i1 = static_cast<signed char>(i1);
                b_i2 = static_cast<signed char>(i3);
                b_i3 = static_cast<signed char>(i2);
                b_i4 = static_cast<signed char>(i4);
              } else {
                b_i1 = static_cast<signed char>(i1);
                b_i2 = static_cast<signed char>(i3);
                b_i3 = static_cast<signed char>(i4);
                b_i4 = static_cast<signed char>(i2);
              }
            } else {
              d1 = x4[i4 - 1];
              if (d <= d1) {
                if (x4[i2 - 1] <= d1) {
                  b_i1 = static_cast<signed char>(i3);
                  b_i2 = static_cast<signed char>(i1);
                  b_i3 = static_cast<signed char>(i2);
                  b_i4 = static_cast<signed char>(i4);
                } else {
                  b_i1 = static_cast<signed char>(i3);
                  b_i2 = static_cast<signed char>(i1);
                  b_i3 = static_cast<signed char>(i4);
                  b_i4 = static_cast<signed char>(i2);
                }
              } else {
                b_i1 = static_cast<signed char>(i3);
                b_i2 = static_cast<signed char>(i4);
                b_i3 = static_cast<signed char>(i1);
                b_i4 = static_cast<signed char>(i2);
              }
            }
            iidx[dim - 3] = idx4[b_i1 - 1];
            iidx[dim - 2] = idx4[b_i2 - 1];
            iidx[dim - 1] = idx4[b_i3 - 1];
            iidx[dim] = idx4[b_i4 - 1];
            vwork[dim - 3] = x4[b_i1 - 1];
            vwork[dim - 2] = x4[b_i2 - 1];
            vwork[dim - 1] = x4[b_i3 - 1];
            vwork[dim] = x4[b_i4 - 1];
            dim = 0;
          }
        }
      }
      i3 = vwork.size(0) - bLen2;
      if (dim > 0) {
        signed char perm[4];
        perm[1] = 0;
        perm[2] = 0;
        perm[3] = 0;
        if (dim == 1) {
          perm[0] = 1;
        } else if (dim == 2) {
          if (x4[0] <= x4[1]) {
            perm[0] = 1;
            perm[1] = 2;
          } else {
            perm[0] = 2;
            perm[1] = 1;
          }
        } else if (x4[0] <= x4[1]) {
          if (x4[1] <= x4[2]) {
            perm[0] = 1;
            perm[1] = 2;
            perm[2] = 3;
          } else if (x4[0] <= x4[2]) {
            perm[0] = 1;
            perm[1] = 3;
            perm[2] = 2;
          } else {
            perm[0] = 3;
            perm[1] = 1;
            perm[2] = 2;
          }
        } else if (x4[0] <= x4[2]) {
          perm[0] = 2;
          perm[1] = 1;
          perm[2] = 3;
        } else if (x4[1] <= x4[2]) {
          perm[0] = 2;
          perm[1] = 3;
          perm[2] = 1;
        } else {
          perm[0] = 3;
          perm[1] = 2;
          perm[2] = 1;
        }
        i = static_cast<unsigned char>(dim);
        for (int k{0}; k < i; k++) {
          iidx_tmp = perm[k] - 1;
          i1 = (i3 - dim) + k;
          iidx[i1] = idx4[iidx_tmp];
          vwork[i1] = x4[iidx_tmp];
        }
      }
      dim = bLen2 >> 1;
      for (int k{0}; k < dim; k++) {
        i1 = i3 + k;
        i2 = iidx[i1];
        iidx_tmp = (n - k) - 1;
        iidx[i1] = iidx[iidx_tmp];
        iidx[iidx_tmp] = i2;
        vwork[i1] = xwork[iidx_tmp];
        vwork[iidx_tmp] = xwork[i1];
      }
      if ((bLen2 & 1) != 0) {
        dim += i3;
        vwork[dim] = xwork[dim];
      }
      nNonNaN = vwork.size(0) - bLen2;
      i1 = 2;
      if (nNonNaN > 1) {
        if (vwork.size(0) >= 256) {
          int nBlocks;
          nBlocks = nNonNaN >> 8;
          if (nBlocks > 0) {
            for (int b{0}; b < nBlocks; b++) {
              double c_xwork[256];
              int c_iwork[256];
              i4 = (b << 8) - 1;
              for (int b_b{0}; b_b < 6; b_b++) {
                n = 1 << (b_b + 2);
                bLen2 = n << 1;
                i = 256 >> (b_b + 3);
                for (int k{0}; k < i; k++) {
                  i2 = (i4 + k * bLen2) + 1;
                  for (i1 = 0; i1 < bLen2; i1++) {
                    dim = i2 + i1;
                    c_iwork[i1] = iidx[dim];
                    c_xwork[i1] = vwork[dim];
                  }
                  i3 = 0;
                  i1 = n;
                  dim = i2 - 1;
                  int exitg1;
                  do {
                    exitg1 = 0;
                    dim++;
                    if (c_xwork[i3] <= c_xwork[i1]) {
                      iidx[dim] = c_iwork[i3];
                      vwork[dim] = c_xwork[i3];
                      if (i3 + 1 < n) {
                        i3++;
                      } else {
                        exitg1 = 1;
                      }
                    } else {
                      iidx[dim] = c_iwork[i1];
                      vwork[dim] = c_xwork[i1];
                      if (i1 + 1 < bLen2) {
                        i1++;
                      } else {
                        dim -= i3;
                        for (i1 = i3 + 1; i1 <= n; i1++) {
                          iidx_tmp = dim + i1;
                          iidx[iidx_tmp] = c_iwork[i1 - 1];
                          vwork[iidx_tmp] = c_xwork[i1 - 1];
                        }
                        exitg1 = 1;
                      }
                    }
                  } while (exitg1 == 0);
                }
              }
            }
            dim = nBlocks << 8;
            i1 = nNonNaN - dim;
            if (i1 > 0) {
              merge_block(iidx, vwork, dim, i1, 2, iwork, xwork);
            }
            i1 = 8;
          }
        }
        dim = iwork.size(0);
        b_iwork.set_size(iwork.size(0));
        for (i = 0; i < dim; i++) {
          b_iwork[i] = iwork[i];
        }
        b_xwork.set_size(xwork.size(0));
        dim = xwork.size(0);
        for (i = 0; i < dim; i++) {
          b_xwork[i] = xwork[i];
        }
        merge_block(iidx, vwork, 0, nNonNaN, i1, b_iwork, b_xwork);
      }
    }
    for (int k{0}; k <= vlen; k++) {
      i = j + k * vstride;
      x[i] = vwork[k];
      idx[i] = iidx[k];
    }
  }
}

//
// Arguments    : double x[4]
// Return Type  : void
//
void sort(double x[4])
{
  double x4[4];
  double xwork[4];
  int i2;
  int i3;
  int ib;
  int nNaNs;
  x4[0] = 0.0;
  xwork[0] = 0.0;
  x4[1] = 0.0;
  xwork[1] = 0.0;
  x4[2] = 0.0;
  xwork[2] = 0.0;
  x4[3] = 0.0;
  xwork[3] = 0.0;
  nNaNs = -3;
  ib = 0;
  if (std::isnan(x[0])) {
    xwork[3] = x[0];
    nNaNs = -2;
  } else {
    ib = 1;
    x4[0] = x[0];
  }
  if (std::isnan(x[1])) {
    xwork[-nNaNs] = x[1];
    nNaNs++;
  } else {
    ib++;
    x4[ib - 1] = x[1];
  }
  if (std::isnan(x[2])) {
    xwork[-nNaNs] = x[2];
    nNaNs++;
  } else {
    ib++;
    x4[ib - 1] = x[2];
  }
  if (std::isnan(x[3])) {
    xwork[-nNaNs] = x[3];
    nNaNs++;
  } else {
    ib++;
    x4[ib - 1] = x[3];
    if (ib == 4) {
      double d;
      double d1;
      int i4;
      signed char b_i2;
      signed char b_i3;
      signed char i;
      signed char i1;
      if (x4[0] <= x4[1]) {
        ib = 1;
        i2 = 2;
      } else {
        ib = 2;
        i2 = 1;
      }
      if (x4[2] <= x4[3]) {
        i3 = 3;
        i4 = 4;
      } else {
        i3 = 4;
        i4 = 3;
      }
      d = x4[ib - 1];
      d1 = x4[i3 - 1];
      if (d <= d1) {
        d = x4[i2 - 1];
        if (d <= d1) {
          i = static_cast<signed char>(ib);
          i1 = static_cast<signed char>(i2);
          b_i2 = static_cast<signed char>(i3);
          b_i3 = static_cast<signed char>(i4);
        } else if (d <= x4[i4 - 1]) {
          i = static_cast<signed char>(ib);
          i1 = static_cast<signed char>(i3);
          b_i2 = static_cast<signed char>(i2);
          b_i3 = static_cast<signed char>(i4);
        } else {
          i = static_cast<signed char>(ib);
          i1 = static_cast<signed char>(i3);
          b_i2 = static_cast<signed char>(i4);
          b_i3 = static_cast<signed char>(i2);
        }
      } else {
        d1 = x4[i4 - 1];
        if (d <= d1) {
          if (x4[i2 - 1] <= d1) {
            i = static_cast<signed char>(i3);
            i1 = static_cast<signed char>(ib);
            b_i2 = static_cast<signed char>(i2);
            b_i3 = static_cast<signed char>(i4);
          } else {
            i = static_cast<signed char>(i3);
            i1 = static_cast<signed char>(ib);
            b_i2 = static_cast<signed char>(i4);
            b_i3 = static_cast<signed char>(i2);
          }
        } else {
          i = static_cast<signed char>(i3);
          i1 = static_cast<signed char>(i4);
          b_i2 = static_cast<signed char>(ib);
          b_i3 = static_cast<signed char>(i2);
        }
      }
      x[-(nNaNs + 3)] = x4[i - 1];
      x[-(nNaNs + 2)] = x4[i1 - 1];
      x[-(nNaNs + 1)] = x4[b_i2 - 1];
      x[-nNaNs] = x4[b_i3 - 1];
      ib = 0;
    }
  }
  if (ib > 0) {
    signed char perm[4];
    perm[1] = 0;
    perm[2] = 0;
    perm[3] = 0;
    if (ib == 1) {
      perm[0] = 1;
    } else if (ib == 2) {
      if (x4[0] <= x4[1]) {
        perm[0] = 1;
        perm[1] = 2;
      } else {
        perm[0] = 2;
        perm[1] = 1;
      }
    } else if (x4[0] <= x4[1]) {
      if (x4[1] <= x4[2]) {
        perm[0] = 1;
        perm[1] = 2;
        perm[2] = 3;
      } else if (x4[0] <= x4[2]) {
        perm[0] = 1;
        perm[1] = 3;
        perm[2] = 2;
      } else {
        perm[0] = 3;
        perm[1] = 1;
        perm[2] = 2;
      }
    } else if (x4[0] <= x4[2]) {
      perm[0] = 2;
      perm[1] = 1;
      perm[2] = 3;
    } else if (x4[1] <= x4[2]) {
      perm[0] = 2;
      perm[1] = 3;
      perm[2] = 1;
    } else {
      perm[0] = 3;
      perm[1] = 2;
      perm[2] = 1;
    }
    for (i3 = 0; i3 < ib; i3++) {
      x[((i3 - nNaNs) - ib) + 1] = x4[perm[i3] - 1];
    }
  }
  ib = ((nNaNs + 3) >> 1) + 1;
  for (i3 = 0; i3 <= ib - 2; i3++) {
    i2 = (i3 - nNaNs) + 1;
    x[i2] = xwork[3 - i3];
    x[3 - i3] = xwork[i2];
  }
  if (((nNaNs + 3) & 1) != 0) {
    i2 = ib - nNaNs;
    x[i2] = xwork[i2];
  }
}

} // namespace internal
} // namespace coder

//
// File trailer for sort.cpp
//
// [EOF]
//
