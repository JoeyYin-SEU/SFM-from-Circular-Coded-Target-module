//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: eml_setop.cpp
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

// Include Files
#include "eml_setop.h"
#include "rt_nonfinite.h"
#include "sort.h"
#include "sortIdx.h"
#include "coder_array.h"
#include "omp.h"
#include <cmath>

// Function Definitions
//
// Arguments    : const ::coder::array<double, 2U> &a
//                const ::coder::array<double, 2U> &b
//                ::coder::array<double, 2U> &c
//                ::coder::array<int, 1U> &ia
//                int *ib_size
// Return Type  : void
//
namespace coder {
void do_vectors(const ::coder::array<double, 2U> &a,
                const ::coder::array<double, 2U> &b,
                ::coder::array<double, 2U> &c, ::coder::array<int, 1U> &ia,
                int *ib_size)
{
  array<int, 2U> aperm;
  array<int, 2U> bperm;
  array<int, 1U> b_ia;
  int b_ialast;
  int iafirst;
  int ialast;
  int iblast;
  int na;
  int nc;
  int nia;
  na = a.size(1);
  c.set_size(1, a.size(1));
  ia.set_size(a.size(1));
  *ib_size = 0;
  internal::sortIdx(a, aperm);
  internal::sortIdx(b, bperm);
  nc = 0;
  nia = -1;
  iafirst = 0;
  ialast = 1;
  iblast = 1;
  while ((ialast <= na) && (iblast <= b.size(1))) {
    double ak;
    double bk;
    b_ialast = ialast;
    ak = a[aperm[ialast - 1] - 1];
    while ((b_ialast < a.size(1)) && (a[aperm[b_ialast] - 1] == ak)) {
      b_ialast++;
    }
    ialast = b_ialast;
    bk = b[bperm[iblast - 1] - 1];
    while ((iblast < b.size(1)) && (b[bperm[iblast] - 1] == bk)) {
      iblast++;
    }
    if (ak == bk) {
      ialast = b_ialast + 1;
      iafirst = b_ialast;
      iblast++;
    } else {
      boolean_T p;
      if (std::isnan(bk)) {
        p = !std::isnan(ak);
      } else if (std::isnan(ak)) {
        p = false;
      } else {
        p = (ak < bk);
      }
      if (p) {
        nc++;
        nia++;
        ia[nia] = aperm[iafirst];
        ialast = b_ialast + 1;
        iafirst = b_ialast;
      } else {
        iblast++;
      }
    }
  }
  while (ialast <= na) {
    b_ialast = ialast;
    while ((b_ialast < a.size(1)) &&
           (a[aperm[b_ialast] - 1] == a[aperm[ialast - 1] - 1])) {
      b_ialast++;
    }
    nc++;
    nia++;
    ia[nia] = aperm[iafirst];
    ialast = b_ialast + 1;
    iafirst = b_ialast;
  }
  if (a.size(1) > 0) {
    if (nia + 1 < 1) {
      na = 0;
    } else {
      na = nia + 1;
    }
    ia.set_size(na);
  }
  internal::sort(ia, b_ia);
  if (static_cast<int>(nia + 1 < 3200)) {
    for (int k{0}; k <= nia; k++) {
      c[k] = a[ia[k] - 1];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k <= nia; k++) {
      c[k] = a[ia[k] - 1];
    }
  }
  if (a.size(1) > 0) {
    if (nc < 1) {
      nc = 0;
    }
    c.set_size(c.size(0), nc);
  }
}

//
// Arguments    : const ::coder::array<double, 2U> &a
//                const ::coder::array<double, 1U> &b
//                ::coder::array<double, 2U> &c
//                ::coder::array<int, 1U> &ia
//                int *ib_size
// Return Type  : void
//
void do_vectors(const ::coder::array<double, 2U> &a,
                const ::coder::array<double, 1U> &b,
                ::coder::array<double, 2U> &c, ::coder::array<int, 1U> &ia,
                int *ib_size)
{
  array<int, 2U> aperm;
  array<int, 1U> b_ia;
  array<int, 1U> bperm;
  array<int, 1U> iwork;
  double ak;
  int i;
  int iafirst;
  int iblast;
  int j;
  int n;
  int na;
  int nc;
  int nia;
  int qEnd;
  unsigned int unnamed_idx_0;
  na = a.size(1);
  c.set_size(1, a.size(1));
  ia.set_size(a.size(1));
  *ib_size = 0;
  internal::sortIdx(a, aperm);
  n = b.size(0) + 1;
  unnamed_idx_0 = static_cast<unsigned int>(b.size(0));
  bperm.set_size(b.size(0));
  iafirst = b.size(0);
  if (static_cast<int>(b.size(0) < 3200)) {
    for (int k{0}; k < iafirst; k++) {
      bperm[k] = 0;
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k < iafirst; k++) {
      bperm[k] = 0;
    }
  }
  if (b.size(0) != 0) {
    int b_k;
    iwork.set_size(static_cast<int>(unnamed_idx_0));
    iafirst = b.size(0) - 1;
    for (b_k = 1; b_k <= iafirst; b_k += 2) {
      if ((b[b_k - 1] <= b[b_k]) || std::isnan(b[b_k])) {
        bperm[b_k - 1] = b_k;
        bperm[b_k] = b_k + 1;
      } else {
        bperm[b_k - 1] = b_k + 1;
        bperm[b_k] = b_k;
      }
    }
    if ((b.size(0) & 1) != 0) {
      bperm[b.size(0) - 1] = b.size(0);
    }
    i = 2;
    while (i < n - 1) {
      iblast = i << 1;
      j = 1;
      for (nc = i + 1; nc < n; nc = qEnd + i) {
        int kEnd;
        int q;
        nia = j;
        q = nc - 1;
        qEnd = j + iblast;
        if (qEnd > n) {
          qEnd = n;
        }
        b_k = 0;
        kEnd = qEnd - j;
        while (b_k + 1 <= kEnd) {
          ak = b[bperm[q] - 1];
          iafirst = bperm[nia - 1];
          if ((b[iafirst - 1] <= ak) || std::isnan(ak)) {
            iwork[b_k] = iafirst;
            nia++;
            if (nia == nc) {
              while (q + 1 < qEnd) {
                b_k++;
                iwork[b_k] = bperm[q];
                q++;
              }
            }
          } else {
            iwork[b_k] = bperm[q];
            q++;
            if (q + 1 == qEnd) {
              while (nia < nc) {
                b_k++;
                iwork[b_k] = bperm[nia - 1];
                nia++;
              }
            }
          }
          b_k++;
        }
        for (b_k = 0; b_k < kEnd; b_k++) {
          bperm[(j + b_k) - 1] = iwork[b_k];
        }
        j = qEnd;
      }
      i = iblast;
    }
  }
  nc = 0;
  nia = -1;
  iafirst = 0;
  i = 1;
  iblast = 1;
  while ((i <= na) && (iblast <= b.size(0))) {
    double bk;
    j = i;
    ak = a[aperm[i - 1] - 1];
    while ((j < a.size(1)) && (a[aperm[j] - 1] == ak)) {
      j++;
    }
    i = j;
    bk = b[bperm[iblast - 1] - 1];
    while ((iblast < b.size(0)) && (b[bperm[iblast] - 1] == bk)) {
      iblast++;
    }
    if (ak == bk) {
      i = j + 1;
      iafirst = j;
      iblast++;
    } else {
      boolean_T p;
      if (std::isnan(bk)) {
        p = !std::isnan(ak);
      } else if (std::isnan(ak)) {
        p = false;
      } else {
        p = (ak < bk);
      }
      if (p) {
        nc++;
        nia++;
        ia[nia] = aperm[iafirst];
        i = j + 1;
        iafirst = j;
      } else {
        iblast++;
      }
    }
  }
  while (i <= na) {
    j = i;
    while ((j < a.size(1)) && (a[aperm[j] - 1] == a[aperm[i - 1] - 1])) {
      j++;
    }
    nc++;
    nia++;
    ia[nia] = aperm[iafirst];
    i = j + 1;
    iafirst = j;
  }
  if (a.size(1) > 0) {
    if (nia + 1 < 1) {
      i = 0;
    } else {
      i = nia + 1;
    }
    ia.set_size(i);
  }
  internal::sort(ia, b_ia);
  if (static_cast<int>(nia + 1 < 3200)) {
    for (int k{0}; k <= nia; k++) {
      c[k] = a[ia[k] - 1];
    }
  } else {
#pragma omp parallel for num_threads(                                          \
    32 > omp_get_max_threads() ? omp_get_max_threads() : 32)

    for (int k = 0; k <= nia; k++) {
      c[k] = a[ia[k] - 1];
    }
  }
  if (a.size(1) > 0) {
    if (nc < 1) {
      nc = 0;
    }
    c.set_size(c.size(0), nc);
  }
}

} // namespace coder

//
// File trailer for eml_setop.cpp
//
// [EOF]
//
