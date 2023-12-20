//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xzsvdc.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 19-Dec-2023 13:39:53
//

// Include Files
#include "xzsvdc.h"
#include "rt_nonfinite.h"
#include "xaxpy.h"
#include "xnrm2.h"
#include "xrotg.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>
#include <cstring>

// Function Definitions
//
// Arguments    : const ::coder::array<double, 2U> &A
//                double *U
//                double S_data[]
//                int *S_size
//                ::coder::array<double, 2U> &V
// Return Type  : void
//
namespace coder {
namespace internal {
namespace reflapack {
void xzsvdc(const ::coder::array<double, 2U> &A, double *U, double S_data[],
            int *S_size, ::coder::array<double, 2U> &V)
{
  array<double, 1U> e;
  double s_data[2];
  double b;
  double c;
  double f;
  double nrm;
  double sm;
  int minnp;
  int nrt;
  int ns;
  int p;
  p = A.size(1);
  if (A.size(1) >= 2) {
    ns = 2;
  } else {
    ns = A.size(1);
  }
  minnp = (A.size(1) >= 1);
  if (A.size(1) >= 2) {
    nrt = 2;
  } else {
    nrt = A.size(1);
  }
  if (static_cast<int>(nrt < 3200)) {
    if (ns - 1 >= 0) {
      std::memset(&s_data[0], 0,
                  static_cast<unsigned int>(ns) * sizeof(double));
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < ns; i++) {
      s_data[i] = 0.0;
    }
  }
  e.set_size(A.size(1));
  ns = A.size(1);
  if (static_cast<int>(A.size(1) < 3200)) {
    for (int i{0}; i < ns; i++) {
      e[i] = 0.0;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < ns; i++) {
      e[i] = 0.0;
    }
  }
  V.set_size(A.size(1), A.size(1));
  ns = A.size(1) * A.size(1);
  if (static_cast<int>(ns < 3200)) {
    for (int i{0}; i < ns; i++) {
      V[i] = 0.0;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int i = 0; i < ns; i++) {
      V[i] = 0.0;
    }
  }
  if (A.size(1) == 0) {
    *U = 1.0;
  } else {
    double rt;
    double snorm;
    int ii;
    int iter;
    int m;
    int mm;
    int q;
    if (A.size(1) >= 2) {
      nrt = A.size(1) - 2;
    } else {
      nrt = 0;
    }
    if (nrt > 1) {
      nrt = 1;
    }
    for (q = 0; q < nrt; q++) {
      for (iter = 2; iter <= p; iter++) {
        e[iter - 1] = A[iter - 1];
      }
      if (nrt >= 1) {
        nrm = blas::xnrm2(p - 1, e, 2);
        if (nrm == 0.0) {
          e[0] = 0.0;
        } else {
          if (e[1] < 0.0) {
            e[0] = -nrm;
          } else {
            e[0] = nrm;
          }
          nrm = e[0];
          if (std::abs(e[0]) >= 1.0020841800044864E-292) {
            nrm = 1.0 / e[0];
            for (int k{2}; k <= p; k++) {
              e[k - 1] = nrm * e[k - 1];
            }
          } else {
            for (int k{2}; k <= p; k++) {
              e[k - 1] = e[k - 1] / nrm;
            }
          }
          e[1] = e[1] + 1.0;
          e[0] = -e[0];
        }
        for (ii = 2; ii <= p; ii++) {
          V[ii - 1] = e[ii - 1];
        }
      }
    }
    m = A.size(1);
    if (m > 2) {
      m = 2;
    }
    s_data[0] = A[0];
    if (m > 1) {
      s_data[1] = 0.0;
    }
    if (nrt + 1 < m) {
      e[0] = A[1];
    }
    e[m - 1] = 0.0;
    *U = 1.0;
    for (q = p; q >= 1; q--) {
      if ((q <= nrt) && (e[0] != 0.0)) {
        for (iter = 2; iter <= p; iter++) {
          ns = p * (iter - 1);
          nrm = 0.0;
          if (p - 1 >= 1) {
            for (int k{0}; k <= p - 2; k++) {
              nrm += V[k + 1] * V[(ns + k) + 1];
            }
          }
          nrm = -(nrm / V[1]);
          blas::xaxpy(p - 1, nrm, V, ns + 2);
        }
      }
      for (ii = 0; ii < p; ii++) {
        V[ii + V.size(0) * (q - 1)] = 0.0;
      }
      V[(q + V.size(0) * (q - 1)) - 1] = 1.0;
    }
    for (q = 0; q < m; q++) {
      nrm = s_data[q];
      if (nrm != 0.0) {
        rt = std::abs(nrm);
        nrm /= rt;
        s_data[q] = rt;
        if (q + 1 < m) {
          e[0] = e[0] / nrm;
        }
        if (q + 1 <= 1) {
          *U *= nrm;
        }
      }
      if ((q + 1 < m) && (e[0] != 0.0)) {
        rt = std::abs(e[0]);
        nrm = rt / e[0];
        e[0] = rt;
        s_data[1] *= nrm;
        ns = p + 1;
        nrt = p + p;
        for (int k{ns}; k <= nrt; k++) {
          V[k - 1] = nrm * V[k - 1];
        }
      }
    }
    mm = m;
    iter = 0;
    snorm = 0.0;
    for (ii = 0; ii < m; ii++) {
      snorm =
          std::fmax(snorm, std::fmax(std::abs(s_data[ii]), std::abs(e[ii])));
    }
    while ((m > 0) && (iter < 75)) {
      boolean_T exitg1;
      q = m - 1;
      ii = m - 1;
      exitg1 = false;
      while (!(exitg1 || (ii == 0))) {
        nrm = std::abs(e[0]);
        if ((nrm <= 2.2204460492503131E-16 *
                        (std::abs(s_data[0]) + std::abs(s_data[1]))) ||
            (nrm <= 1.0020841800044864E-292) ||
            ((iter > 20) && (nrm <= 2.2204460492503131E-16 * snorm))) {
          e[0] = 0.0;
          exitg1 = true;
        } else {
          ii = 0;
        }
      }
      if (ii == m - 1) {
        ns = 4;
      } else {
        nrt = m;
        ns = m;
        exitg1 = false;
        while ((!exitg1) && (ns >= ii)) {
          nrt = ns;
          if (ns == ii) {
            exitg1 = true;
          } else {
            nrm = 0.0;
            if (ns < m) {
              nrm = std::abs(e[0]);
            }
            if (ns > ii + 1) {
              nrm += std::abs(e[0]);
            }
            rt = std::abs(s_data[ns - 1]);
            if ((rt <= 2.2204460492503131E-16 * nrm) ||
                (rt <= 1.0020841800044864E-292)) {
              s_data[ns - 1] = 0.0;
              exitg1 = true;
            } else {
              ns--;
            }
          }
        }
        if (nrt == ii) {
          ns = 3;
        } else if (nrt == m) {
          ns = 1;
        } else {
          ns = 2;
          ii = nrt;
        }
      }
      switch (ns) {
      case 1:
        f = e[0];
        e[0] = 0.0;
        for (int k{q}; k >= ii + 1; k--) {
          blas::xrotg(&s_data[0], &f, &c, &nrm);
          if (p >= 1) {
            ns = p * (m - 1);
            for (int b_k{0}; b_k < p; b_k++) {
              nrt = ns + b_k;
              rt = c * V[b_k] + nrm * V[nrt];
              V[nrt] = c * V[nrt] - nrm * V[b_k];
              V[b_k] = rt;
            }
          }
        }
        break;
      case 2:
        f = e[ii - 1];
        e[ii - 1] = 0.0;
        for (int k{ii + 1}; k <= m; k++) {
          blas::xrotg(&s_data[k - 1], &f, &b, &sm);
          nrm = e[k - 1];
          f = -sm * nrm;
          e[k - 1] = nrm * b;
          *U = b * *U + sm * *U;
        }
        break;
      case 3: {
        double scale;
        double sqds;
        nrm = s_data[m - 1];
        scale = std::fmax(
            std::fmax(std::fmax(std::fmax(std::abs(nrm), std::abs(s_data[0])),
                                std::abs(e[0])),
                      std::abs(s_data[ii])),
            std::abs(e[ii]));
        sm = nrm / scale;
        nrm = s_data[0] / scale;
        rt = e[0] / scale;
        sqds = s_data[ii] / scale;
        b = ((nrm + sm) * (nrm - sm) + rt * rt) / 2.0;
        c = sm * rt;
        c *= c;
        if ((b != 0.0) || (c != 0.0)) {
          nrm = std::sqrt(b * b + c);
          if (b < 0.0) {
            nrm = -nrm;
          }
          nrm = c / (b + nrm);
        } else {
          nrm = 0.0;
        }
        f = (sqds + sm) * (sqds - sm) + nrm;
        nrm = sqds * (e[ii] / scale);
        for (int k{ii + 1}; k < 2; k++) {
          blas::xrotg(&f, &nrm, &b, &sm);
          f = b * s_data[0] + sm * e[0];
          e[0] = b * e[0] - sm * s_data[0];
          nrm = sm * s_data[1];
          s_data[1] *= b;
          if (p >= 1) {
            for (int b_k{0}; b_k < p; b_k++) {
              rt = b * V[b_k] + sm * V[p + b_k];
              V[p + b_k] = b * V[p + b_k] - sm * V[b_k];
              V[b_k] = rt;
            }
          }
          s_data[0] = f;
          blas::xrotg(&s_data[0], &nrm, &b, &sm);
          f = b * e[0] + sm * s_data[1];
          s_data[1] = -sm * e[0] + b * s_data[1];
          nrm = sm * e[1];
          e[1] = e[1] * b;
        }
        e[0] = f;
        iter++;
      } break;
      default:
        if (s_data[ii] < 0.0) {
          s_data[ii] = -s_data[ii];
          ns = p * ii;
          nrt = ns + p;
          for (int k{ns + 1}; k <= nrt; k++) {
            V[k - 1] = -V[k - 1];
          }
        }
        while ((ii + 1 < mm) && (s_data[0] < s_data[1])) {
          rt = s_data[0];
          s_data[0] = s_data[1];
          s_data[1] = rt;
          if (p > 1) {
            for (int k{0}; k < p; k++) {
              rt = V[k];
              V[k] = V[p + k];
              V[p + k] = rt;
            }
          }
          ii = 1;
        }
        iter = 0;
        m--;
        break;
      }
    }
  }
  *S_size = minnp;
  if (minnp - 1 >= 0) {
    S_data[0] = s_data[0];
  }
}

} // namespace reflapack
} // namespace internal
} // namespace coder

//
// File trailer for xzsvdc.cpp
//
// [EOF]
//
