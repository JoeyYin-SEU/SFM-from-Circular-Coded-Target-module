//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: Checkerboard.h
//
// MATLAB Coder version            : 5.5
// C/C++ source code generated on  : 27-Nov-2023 10:57:33
//

#ifndef CHECKERBOARD_H
#define CHECKERBOARD_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Type Definitions
namespace coder {
namespace vision {
namespace internal {
namespace calibration {
namespace checkerboard {
class Checkerboard {
public:
  void initialize(double seedIdx, const ::coder::array<float, 2U> &points,
                  const float v1[2], const float v2[2]);
  double findNeighbor(const ::coder::array<float, 2U> &pointVectors,
                      const ::coder::array<float, 1U> &euclideanDists,
                      const ::coder::array<float, 2U> &v) const;
  boolean_T expandBoardOnce();
  void findClosestIndices(const ::coder::array<double, 2U> &predictedPoints,
                          ::coder::array<double, 2U> &indices) const;
  boolean_T b_expandBoardOnce();
  static void arrayFind(const ::coder::array<boolean_T, 2U> &arr,
                        ::coder::array<double, 2U> &matchedIdx);

private:
  double findNeighbor(const ::coder::array<float, 2U> &pointVectors,
                      const ::coder::array<float, 1U> &euclideanDists,
                      const float v[2]) const;
  float computeInitialEnergy() const;
  void fitPolynomialIndices(::coder::array<double, 2U> &newIndices) const;
  void findIndependentVar(double coordsToUse[2]) const;
  void findClosestOnCurve(const double predictedPoint[2], double radius,
                          const double curve_data[], const int curve_size[2],
                          const double coordsToUse[2],
                          const ::coder::array<double, 2U> &removedIdx,
                          ::coder::array<double, 2U> &idx) const;
  void expandBoardUp(const ::coder::array<double, 2U> &indices,
                     ::coder::array<double, 2U> &newBoard,
                     ::coder::array<double, 3U> &newBoardCoords) const;
  float computeNewEnergyVertical(float oldEnergy) const;
  void findIndependentVar(const ::coder::array<double, 2U> &idx,
                          double coordsToUse[2]) const;
  void expandBoardDown(const ::coder::array<double, 2U> &indices,
                       ::coder::array<double, 2U> &newBoard,
                       ::coder::array<double, 3U> &newBoardCoords) const;
  float computeNewEnergyVertical(const ::coder::array<double, 2U> &idx,
                                 float oldEnergy) const;
  void b_fitPolynomialIndices(::coder::array<double, 2U> &newIndices) const;
  void b_findIndependentVar(double coordsToUse[2]) const;
  void expandBoardLeft(const ::coder::array<double, 2U> &indices,
                       ::coder::array<double, 2U> &newBoard,
                       ::coder::array<double, 3U> &newBoardCoords) const;
  float computeNewEnergyHorizontal(float oldEnergy) const;
  void fitPolynomialIndices(const ::coder::array<double, 2U> &idx,
                            ::coder::array<double, 2U> &newIndices) const;
  void b_findIndependentVar(const ::coder::array<double, 2U> &idx,
                            double coordsToUse[2]) const;
  void expandBoardRight(const ::coder::array<double, 2U> &indices,
                        ::coder::array<double, 2U> &newBoard,
                        ::coder::array<double, 3U> &newBoardCoords) const;
  float computeNewEnergyHorizontal(const ::coder::array<double, 2U> &idx,
                                   float oldEnergy) const;
  void undoLastExpansion();
  void c_fitPolynomialIndices(::coder::array<double, 2U> &newIndices) const;
  void b_findClosestIndices(const ::coder::array<double, 2U> &predictedPoints,
                            ::coder::array<double, 2U> &indices) const;
  void d_fitPolynomialIndices(::coder::array<double, 2U> &newIndices) const;
  void b_fitPolynomialIndices(const ::coder::array<double, 2U> &idx,
                              ::coder::array<double, 2U> &newIndices) const;
  void findSearchParams(const ::coder::array<double, 2U> &idx,
                        const ::coder::array<double, 1U> &validIdx,
                        double currIdx, const double coordsToUse[2],
                        double *coordDist, double *moveMultiplier,
                        double *firstValidIdx) const;
  void findSearchParams(const ::coder::array<double, 2U> &idx,
                        const ::coder::array<double, 2U> &validIdx,
                        double currIdx, const double coordsToUse[2],
                        double *coordDist, double *moveMultiplier,
                        double *firstValidIdx) const;

public:
  boolean_T isValid;
  float Energy;
  array<double, 3U> BoardCoords;
  array<double, 2U> BoardIdx;
  array<float, 2U> Points;
  boolean_T IsDirectionBad[4];
  boolean_T IsDistortionHigh;

private:
  double LastExpandDirection;
  float PreviousEnergy;
};

} // namespace checkerboard
} // namespace calibration
} // namespace internal
} // namespace vision
} // namespace coder

#endif
//
// File trailer for Checkerboard.h
//
// [EOF]
//
