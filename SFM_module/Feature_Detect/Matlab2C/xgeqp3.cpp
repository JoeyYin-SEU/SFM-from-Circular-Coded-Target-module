//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xgeqp3.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "xgeqp3.h"
#include "get_chessborad_pixel_rtwutil.h"
#include "rt_nonfinite.h"
#include "xnrm2.h"
#include "coder_array.h"
#include <cmath>
#include <cstring>

// Function Definitions
//
// Arguments    : ::coder::array<double, 2U> &A
//                double tau_data[]
//                int *tau_size
//                int jpvt_data[]
//                int jpvt_size[2]
// Return Type  : void
//
namespace coder {
namespace internal {
namespace lapack {
void xgeqp3(::coder::array<double, 2U> &A, double tau_data[], int *tau_size,
            int jpvt_data[], int jpvt_size[2])
{
  array<double, 2U> x;
  double vn1_data[5];
  double vn2_data[5];
  double work_data[5];
  int itemp;
  int ix;
  int lastc;
  int m;
  int n;
  boolean_T guard1{false};
  m = A.size(0);
  n = A.size(1);
  itemp = A.size(0);
  *tau_size = A.size(1);
  if (itemp <= *tau_size) {
    *tau_size = itemp;
  }
  if (*tau_size - 1 >= 0) {
    std::memset(&tau_data[0], 0,
                static_cast<unsigned int>(*tau_size) * sizeof(double));
  }
  guard1 = false;
  if (A.size(0) == 0) {
    guard1 = true;
  } else {
    int u1;
    itemp = A.size(0);
    u1 = A.size(1);
    if (itemp <= u1) {
      u1 = itemp;
    }
    if (u1 < 1) {
      guard1 = true;
    } else {
      double absxk;
      int ma;
      jpvt_size[0] = 1;
      jpvt_size[1] = A.size(1);
      lastc = A.size(1);
      if (lastc - 1 >= 0) {
        std::memset(&jpvt_data[0], 0,
                    static_cast<unsigned int>(lastc) * sizeof(int));
      }
      for (lastc = 0; lastc < n; lastc++) {
        jpvt_data[lastc] = lastc + 1;
      }
      ma = A.size(0);
      itemp = A.size(0);
      u1 = A.size(1);
      if (itemp <= u1) {
        u1 = itemp;
      }
      lastc = A.size(1);
      if (lastc - 1 >= 0) {
        std::memset(&work_data[0], 0,
                    static_cast<unsigned int>(lastc) * sizeof(double));
      }
      lastc = A.size(1);
      if (lastc - 1 >= 0) {
        std::memset(&vn1_data[0], 0,
                    static_cast<unsigned int>(lastc) * sizeof(double));
      }
      lastc = A.size(1);
      if (lastc - 1 >= 0) {
        std::memset(&vn2_data[0], 0,
                    static_cast<unsigned int>(lastc) * sizeof(double));
      }
      for (ix = 0; ix < n; ix++) {
        absxk = blas::xnrm2(m, A, ix * ma + 1);
        vn1_data[ix] = absxk;
        vn2_data[ix] = absxk;
      }
      for (int i{0}; i < u1; i++) {
        double s;
        double smax;
        double t;
        int b_i;
        int i1;
        int ii;
        int ii_tmp;
        int ip1;
        int mmi;
        int nmi;
        int pvt;
        ip1 = i + 2;
        ii_tmp = i * ma;
        ii = ii_tmp + i;
        nmi = n - i;
        mmi = (m - i) - 1;
        if (nmi < 1) {
          itemp = -1;
        } else {
          itemp = 0;
          if (nmi > 1) {
            smax = std::abs(vn1_data[i]);
            for (lastc = 2; lastc <= nmi; lastc++) {
              s = std::abs(vn1_data[(i + lastc) - 1]);
              if (s > smax) {
                itemp = lastc - 1;
                smax = s;
              }
            }
          }
        }
        pvt = i + itemp;
        if (pvt != i) {
          ix = pvt * ma;
          x.set_size(A.size(0), A.size(1));
          lastc = A.size(0) * A.size(1);
          for (b_i = 0; b_i < lastc; b_i++) {
            x[b_i] = A[b_i];
          }
          for (lastc = 0; lastc < m; lastc++) {
            itemp = ix + lastc;
            smax = x[itemp];
            b_i = ii_tmp + lastc;
            x[itemp] = x[b_i];
            x[b_i] = smax;
          }
          A.set_size(x.size(0), x.size(1));
          lastc = x.size(1);
          for (b_i = 0; b_i < lastc; b_i++) {
            itemp = x.size(0);
            for (i1 = 0; i1 < itemp; i1++) {
              A[i1 + A.size(0) * b_i] = x[i1 + x.size(0) * b_i];
            }
          }
          itemp = jpvt_data[pvt];
          jpvt_data[pvt] = jpvt_data[i];
          jpvt_data[i] = itemp;
          vn1_data[pvt] = vn1_data[i];
          vn2_data[pvt] = vn2_data[i];
        }
        if (i + 1 < m) {
          t = A[ii];
          pvt = ii + 2;
          tau_data[i] = 0.0;
          if (mmi + 1 > 0) {
            smax = blas::xnrm2(mmi, A, ii + 2);
            if (smax != 0.0) {
              s = rt_hypotd_snf(A[ii], smax);
              if (A[ii] >= 0.0) {
                s = -s;
              }
              if (std::abs(s) < 1.0020841800044864E-292) {
                ix = 0;
                b_i = ii + mmi;
                do {
                  ix++;
                  x.set_size(A.size(0), A.size(1));
                  lastc = A.size(0) * A.size(1);
                  for (i1 = 0; i1 < lastc; i1++) {
                    x[i1] = A[i1];
                  }
                  for (lastc = pvt; lastc <= b_i + 1; lastc++) {
                    x[lastc - 1] = 9.9792015476736E+291 * x[lastc - 1];
                  }
                  A.set_size(x.size(0), x.size(1));
                  lastc = x.size(1);
                  for (i1 = 0; i1 < lastc; i1++) {
                    itemp = x.size(0);
                    for (ii_tmp = 0; ii_tmp < itemp; ii_tmp++) {
                      A[ii_tmp + A.size(0) * i1] = x[ii_tmp + x.size(0) * i1];
                    }
                  }
                  s *= 9.9792015476736E+291;
                  t *= 9.9792015476736E+291;
                } while ((std::abs(s) < 1.0020841800044864E-292) && (ix < 20));
                s = rt_hypotd_snf(t, blas::xnrm2(mmi, x, ii + 2));
                if (t >= 0.0) {
                  s = -s;
                }
                tau_data[i] = (s - t) / s;
                smax = 1.0 / (t - s);
                b_i = ii + mmi;
                for (lastc = pvt; lastc <= b_i + 1; lastc++) {
                  x[lastc - 1] = smax * x[lastc - 1];
                }
                A.set_size(x.size(0), x.size(1));
                lastc = x.size(1);
                for (b_i = 0; b_i < lastc; b_i++) {
                  itemp = x.size(0);
                  for (i1 = 0; i1 < itemp; i1++) {
                    A[i1 + A.size(0) * b_i] = x[i1 + x.size(0) * b_i];
                  }
                }
                for (lastc = 0; lastc < ix; lastc++) {
                  s *= 1.0020841800044864E-292;
                }
                t = s;
              } else {
                tau_data[i] = (s - A[ii]) / s;
                smax = 1.0 / (A[ii] - s);
                x.set_size(A.size(0), A.size(1));
                lastc = A.size(0) * A.size(1);
                for (b_i = 0; b_i < lastc; b_i++) {
                  x[b_i] = A[b_i];
                }
                b_i = ii + mmi;
                for (lastc = pvt; lastc <= b_i + 1; lastc++) {
                  x[lastc - 1] = smax * x[lastc - 1];
                }
                A.set_size(x.size(0), x.size(1));
                lastc = x.size(1);
                for (b_i = 0; b_i < lastc; b_i++) {
                  itemp = x.size(0);
                  for (i1 = 0; i1 < itemp; i1++) {
                    A[i1 + A.size(0) * b_i] = x[i1 + x.size(0) * b_i];
                  }
                }
                t = s;
              }
            }
          }
          A[ii] = t;
        } else {
          tau_data[i] = 0.0;
        }
        if (i + 1 < n) {
          int jA;
          t = A[ii];
          A[ii] = 1.0;
          jA = (ii + ma) + 1;
          if (tau_data[i] != 0.0) {
            boolean_T exitg2;
            pvt = mmi;
            itemp = ii + mmi;
            while ((pvt + 1 > 0) && (A[itemp] == 0.0)) {
              pvt--;
              itemp--;
            }
            lastc = nmi - 2;
            exitg2 = false;
            while ((!exitg2) && (lastc + 1 > 0)) {
              int exitg1;
              itemp = jA + lastc * ma;
              ii_tmp = itemp;
              do {
                exitg1 = 0;
                if (ii_tmp <= itemp + pvt) {
                  if (A[ii_tmp - 1] != 0.0) {
                    exitg1 = 1;
                  } else {
                    ii_tmp++;
                  }
                } else {
                  lastc--;
                  exitg1 = 2;
                }
              } while (exitg1 == 0);
              if (exitg1 == 1) {
                exitg2 = true;
              }
            }
          } else {
            pvt = -1;
            lastc = -1;
          }
          if (pvt + 1 > 0) {
            if (lastc + 1 != 0) {
              if (lastc >= 0) {
                std::memset(&work_data[0], 0,
                            static_cast<unsigned int>(lastc + 1) *
                                sizeof(double));
              }
              itemp = 0;
              b_i = jA + ma * lastc;
              for (ix = jA; ma < 0 ? ix >= b_i : ix <= b_i; ix += ma) {
                smax = 0.0;
                i1 = ix + pvt;
                for (ii_tmp = ix; ii_tmp <= i1; ii_tmp++) {
                  smax += A[ii_tmp - 1] * A[(ii + ii_tmp) - ix];
                }
                work_data[itemp] += smax;
                itemp++;
              }
            }
            if (!(-tau_data[i] == 0.0)) {
              for (ix = 0; ix <= lastc; ix++) {
                absxk = work_data[ix];
                if (absxk != 0.0) {
                  smax = absxk * -tau_data[i];
                  b_i = pvt + jA;
                  for (itemp = jA; itemp <= b_i; itemp++) {
                    A[itemp - 1] = A[itemp - 1] + A[(ii + itemp) - jA] * smax;
                  }
                }
                jA += ma;
              }
            }
          }
          A[ii] = t;
        }
        for (ix = ip1; ix <= n; ix++) {
          itemp = (i + (ix - 1) * ma) + 1;
          absxk = vn1_data[ix - 1];
          if (absxk != 0.0) {
            smax = std::abs(A[itemp - 1]) / absxk;
            smax = 1.0 - smax * smax;
            if (smax < 0.0) {
              smax = 0.0;
            }
            s = absxk / vn2_data[ix - 1];
            s = smax * (s * s);
            if (s <= 1.4901161193847656E-8) {
              if (i + 1 < m) {
                pvt = itemp + 1;
                smax = 0.0;
                if (mmi >= 1) {
                  if (mmi == 1) {
                    smax = std::abs(A[itemp]);
                  } else {
                    s = 3.3121686421112381E-170;
                    itemp += mmi;
                    for (lastc = pvt; lastc <= itemp; lastc++) {
                      absxk = std::abs(A[lastc - 1]);
                      if (absxk > s) {
                        t = s / absxk;
                        smax = smax * t * t + 1.0;
                        s = absxk;
                      } else {
                        t = absxk / s;
                        smax += t * t;
                      }
                    }
                    smax = s * std::sqrt(smax);
                  }
                }
                vn1_data[ix - 1] = smax;
                vn2_data[ix - 1] = smax;
              } else {
                vn1_data[ix - 1] = 0.0;
                vn2_data[ix - 1] = 0.0;
              }
            } else {
              vn1_data[ix - 1] = absxk * std::sqrt(smax);
            }
          }
        }
      }
    }
  }
  if (guard1) {
    jpvt_size[0] = 1;
    jpvt_size[1] = A.size(1);
    lastc = A.size(1);
    if (lastc - 1 >= 0) {
      std::memset(&jpvt_data[0], 0,
                  static_cast<unsigned int>(lastc) * sizeof(int));
    }
    for (ix = 0; ix < n; ix++) {
      jpvt_data[ix] = ix + 1;
    }
  }
}

} // namespace lapack
} // namespace internal
} // namespace coder

//
// File trailer for xgeqp3.cpp
//
// [EOF]
//
