//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: svd1.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

// Include Files
#include "svd1.h"
#include "rt_nonfinite.h"
#include "xaxpy.h"
#include "xdotc.h"
#include "xnrm2.h"
#include "xrot.h"
#include "xrotg.h"
#include "xswap.h"
#include <algorithm>
#include <cmath>
#include <cstring>

// Function Definitions
//
// Arguments    : const double A[841]
//                double U[841]
//                double s[29]
//                double V[841]
// Return Type  : void
//
namespace coder {
namespace internal {
void b_svd(const double A[841], double U[841], double s[29], double V[841])
{
  double b_A[841];
  double e[29];
  double work[29];
  double nrm;
  double rt;
  double sm;
  double snorm;
  double sqds;
  int ii;
  int jj;
  int m;
  int qjj;
  int qp1;
  int qp1jj;
  int qq;
  std::copy(&A[0], &A[841], &b_A[0]);
  std::memset(&s[0], 0, 29U * sizeof(double));
  std::memset(&e[0], 0, 29U * sizeof(double));
  std::memset(&work[0], 0, 29U * sizeof(double));
  std::memset(&U[0], 0, 841U * sizeof(double));
  std::memset(&V[0], 0, 841U * sizeof(double));
  for (int q{0}; q < 28; q++) {
    boolean_T apply_transform;
    qp1 = q + 2;
    qq = (q + 29 * q) + 1;
    apply_transform = false;
    nrm = blas::xnrm2(29 - q, b_A, qq);
    if (nrm > 0.0) {
      apply_transform = true;
      if (b_A[qq - 1] < 0.0) {
        nrm = -nrm;
      }
      s[q] = nrm;
      if (std::abs(nrm) >= 1.0020841800044864E-292) {
        nrm = 1.0 / nrm;
        qp1jj = (qq - q) + 28;
        for (int k{qq}; k <= qp1jj; k++) {
          b_A[k - 1] *= nrm;
        }
      } else {
        qp1jj = (qq - q) + 28;
        for (int k{qq}; k <= qp1jj; k++) {
          b_A[k - 1] /= s[q];
        }
      }
      b_A[qq - 1]++;
      s[q] = -s[q];
    } else {
      s[q] = 0.0;
    }
    for (jj = qp1; jj < 30; jj++) {
      qjj = q + 29 * (jj - 1);
      if (apply_transform) {
        blas::xaxpy(
            29 - q,
            -(blas::xdotc(29 - q, b_A, qq, b_A, qjj + 1) / b_A[q + 29 * q]), qq,
            b_A, qjj + 1);
      }
      e[jj - 1] = b_A[qjj];
    }
    for (ii = q + 1; ii < 30; ii++) {
      qjj = (ii + 29 * q) - 1;
      U[qjj] = b_A[qjj];
    }
    if (q + 1 <= 27) {
      nrm = blas::b_xnrm2(28 - q, e, q + 2);
      if (nrm == 0.0) {
        e[q] = 0.0;
      } else {
        if (e[q + 1] < 0.0) {
          e[q] = -nrm;
        } else {
          e[q] = nrm;
        }
        nrm = e[q];
        if (std::abs(e[q]) >= 1.0020841800044864E-292) {
          nrm = 1.0 / e[q];
          for (int k{qp1}; k < 30; k++) {
            e[k - 1] *= nrm;
          }
        } else {
          for (int k{qp1}; k < 30; k++) {
            e[k - 1] /= nrm;
          }
        }
        e[q + 1]++;
        e[q] = -e[q];
        for (ii = qp1; ii < 30; ii++) {
          work[ii - 1] = 0.0;
        }
        for (jj = qp1; jj < 30; jj++) {
          blas::xaxpy(28 - q, e[jj - 1], b_A, (q + 29 * (jj - 1)) + 2, work,
                      q + 2);
        }
        for (jj = qp1; jj < 30; jj++) {
          blas::b_xaxpy(28 - q, -e[jj - 1] / e[q + 1], work, q + 2, b_A,
                        (q + 29 * (jj - 1)) + 2);
        }
      }
      for (ii = qp1; ii < 30; ii++) {
        V[(ii + 29 * q) - 1] = e[ii - 1];
      }
    }
  }
  m = 27;
  s[28] = b_A[840];
  e[27] = b_A[839];
  e[28] = 0.0;
  std::memset(&U[812], 0, 29U * sizeof(double));
  U[840] = 1.0;
  for (int q{27}; q >= 0; q--) {
    qp1 = q + 2;
    qq = q + 29 * q;
    if (s[q] != 0.0) {
      for (jj = qp1; jj < 30; jj++) {
        qjj = (q + 29 * (jj - 1)) + 1;
        blas::xaxpy(29 - q, -(blas::xdotc(29 - q, U, qq + 1, U, qjj) / U[qq]),
                    qq + 1, U, qjj);
      }
      for (ii = q + 1; ii < 30; ii++) {
        qjj = (ii + 29 * q) - 1;
        U[qjj] = -U[qjj];
      }
      U[qq]++;
      for (ii = 0; ii < q; ii++) {
        U[ii + 29 * q] = 0.0;
      }
    } else {
      std::memset(&U[q * 29], 0, 29U * sizeof(double));
      U[qq] = 1.0;
    }
  }
  for (int q{28}; q >= 0; q--) {
    if ((q + 1 <= 27) && (e[q] != 0.0)) {
      qp1 = q + 2;
      qjj = (q + 29 * q) + 2;
      for (jj = qp1; jj < 30; jj++) {
        qp1jj = (q + 29 * (jj - 1)) + 2;
        blas::xaxpy(28 - q,
                    -(blas::xdotc(28 - q, V, qjj, V, qp1jj) / V[qjj - 1]), qjj,
                    V, qp1jj);
      }
    }
    std::memset(&V[q * 29], 0, 29U * sizeof(double));
    V[q + 29 * q] = 1.0;
  }
  qq = 0;
  snorm = 0.0;
  for (int q{0}; q < 29; q++) {
    nrm = s[q];
    if (nrm != 0.0) {
      rt = std::abs(nrm);
      nrm /= rt;
      s[q] = rt;
      if (q + 1 < 29) {
        e[q] /= nrm;
      }
      qjj = 29 * q;
      qp1jj = qjj + 29;
      for (int k{qjj + 1}; k <= qp1jj; k++) {
        U[k - 1] *= nrm;
      }
    }
    if (q + 1 < 29) {
      nrm = e[q];
      if (nrm != 0.0) {
        rt = std::abs(nrm);
        nrm = rt / nrm;
        e[q] = rt;
        s[q + 1] *= nrm;
        qjj = 29 * (q + 1);
        qp1jj = qjj + 29;
        for (int k{qjj + 1}; k <= qp1jj; k++) {
          V[k - 1] *= nrm;
        }
      }
    }
    snorm = std::fmax(snorm, std::fmax(std::abs(s[q]), std::abs(e[q])));
  }
  while ((m + 2 > 0) && (qq < 75)) {
    boolean_T exitg1;
    jj = m + 1;
    ii = m + 1;
    exitg1 = false;
    while (!(exitg1 || (ii == 0))) {
      nrm = std::abs(e[ii - 1]);
      if ((nrm <=
           2.2204460492503131E-16 * (std::abs(s[ii - 1]) + std::abs(s[ii]))) ||
          (nrm <= 1.0020841800044864E-292) ||
          ((qq > 20) && (nrm <= 2.2204460492503131E-16 * snorm))) {
        e[ii - 1] = 0.0;
        exitg1 = true;
      } else {
        ii--;
      }
    }
    if (ii == m + 1) {
      qjj = 4;
    } else {
      qp1jj = m + 2;
      qjj = m + 2;
      exitg1 = false;
      while ((!exitg1) && (qjj >= ii)) {
        qp1jj = qjj;
        if (qjj == ii) {
          exitg1 = true;
        } else {
          nrm = 0.0;
          if (qjj < m + 2) {
            nrm = std::abs(e[qjj - 1]);
          }
          if (qjj > ii + 1) {
            nrm += std::abs(e[qjj - 2]);
          }
          rt = std::abs(s[qjj - 1]);
          if ((rt <= 2.2204460492503131E-16 * nrm) ||
              (rt <= 1.0020841800044864E-292)) {
            s[qjj - 1] = 0.0;
            exitg1 = true;
          } else {
            qjj--;
          }
        }
      }
      if (qp1jj == ii) {
        qjj = 3;
      } else if (qp1jj == m + 2) {
        qjj = 1;
      } else {
        qjj = 2;
        ii = qp1jj;
      }
    }
    switch (qjj) {
    case 1: {
      rt = e[m];
      e[m] = 0.0;
      for (int k{jj}; k >= ii + 1; k--) {
        blas::xrotg(&s[k - 1], &rt, &sm, &sqds);
        if (k > ii + 1) {
          double b;
          b = e[k - 2];
          rt = -sqds * b;
          e[k - 2] = b * sm;
        }
        blas::xrot(V, 29 * (k - 1) + 1, 29 * (m + 1) + 1, sm, sqds);
      }
    } break;
    case 2: {
      rt = e[ii - 1];
      e[ii - 1] = 0.0;
      for (int k{ii + 1}; k <= m + 2; k++) {
        double b;
        blas::xrotg(&s[k - 1], &rt, &sm, &sqds);
        b = e[k - 1];
        rt = -sqds * b;
        e[k - 1] = b * sm;
        blas::xrot(U, 29 * (k - 1) + 1, 29 * (ii - 1) + 1, sm, sqds);
      }
    } break;
    case 3: {
      double b;
      double scale;
      nrm = s[m + 1];
      scale = std::fmax(
          std::fmax(std::fmax(std::fmax(std::abs(nrm), std::abs(s[m])),
                              std::abs(e[m])),
                    std::abs(s[ii])),
          std::abs(e[ii]));
      sm = nrm / scale;
      nrm = s[m] / scale;
      rt = e[m] / scale;
      sqds = s[ii] / scale;
      b = ((nrm + sm) * (nrm - sm) + rt * rt) / 2.0;
      nrm = sm * rt;
      nrm *= nrm;
      if ((b != 0.0) || (nrm != 0.0)) {
        rt = std::sqrt(b * b + nrm);
        if (b < 0.0) {
          rt = -rt;
        }
        rt = nrm / (b + rt);
      } else {
        rt = 0.0;
      }
      rt += (sqds + sm) * (sqds - sm);
      nrm = sqds * (e[ii] / scale);
      for (int k{ii + 1}; k <= jj; k++) {
        blas::xrotg(&rt, &nrm, &sm, &sqds);
        if (k > ii + 1) {
          e[k - 2] = rt;
        }
        nrm = e[k - 1];
        b = s[k - 1];
        e[k - 1] = sm * nrm - sqds * b;
        rt = sqds * s[k];
        s[k] *= sm;
        blas::xrot(V, 29 * (k - 1) + 1, 29 * k + 1, sm, sqds);
        s[k - 1] = sm * b + sqds * nrm;
        blas::xrotg(&s[k - 1], &rt, &sm, &sqds);
        rt = sm * e[k - 1] + sqds * s[k];
        s[k] = -sqds * e[k - 1] + sm * s[k];
        nrm = sqds * e[k];
        e[k] *= sm;
        blas::xrot(U, 29 * (k - 1) + 1, 29 * k + 1, sm, sqds);
      }
      e[m] = rt;
      qq++;
    } break;
    default:
      if (s[ii] < 0.0) {
        s[ii] = -s[ii];
        qjj = 29 * ii;
        qp1jj = qjj + 29;
        for (int k{qjj + 1}; k <= qp1jj; k++) {
          V[k - 1] = -V[k - 1];
        }
      }
      qp1 = ii + 1;
      while ((ii + 1 < 29) && (s[ii] < s[qp1])) {
        rt = s[ii];
        s[ii] = s[qp1];
        s[qp1] = rt;
        blas::xswap(V, 29 * ii + 1, 29 * (ii + 1) + 1);
        blas::xswap(U, 29 * ii + 1, 29 * (ii + 1) + 1);
        ii = qp1;
        qp1++;
      }
      qq = 0;
      m--;
      break;
    }
  }
}

} // namespace internal
} // namespace coder

//
// File trailer for svd1.cpp
//
// [EOF]
//
