//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: interp1.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 05-Aug-2021 13:41:18
//

// Include Files
#include "interp1.h"
#include "rt_nonfinite.h"
#include <algorithm>
#include <cmath>

// Function Definitions
//
// Arguments    : const double varargin_1[6]
//                const double varargin_2[6]
//                double varargin_3
// Return Type  : double
//
namespace coder {
double b_interp1(const double varargin_1[6], const double varargin_2[6],
                 double varargin_3)
{
  double x[6];
  double y[6];
  double Vq;
  int low_i;
  for (low_i = 0; low_i < 6; low_i++) {
    y[low_i] = varargin_2[low_i];
    x[low_i] = varargin_1[low_i];
  }
  low_i = 0;
  int exitg1;
  do {
    exitg1 = 0;
    if (low_i < 6) {
      if (std::isnan(varargin_1[low_i])) {
        exitg1 = 1;
      } else {
        low_i++;
      }
    } else {
      double xtmp;
      if (varargin_1[1] < varargin_1[0]) {
        xtmp = x[0];
        x[0] = x[5];
        x[5] = xtmp;
        xtmp = y[0];
        y[0] = y[5];
        y[5] = xtmp;
        xtmp = x[1];
        x[1] = x[4];
        x[4] = xtmp;
        xtmp = y[1];
        y[1] = y[4];
        y[4] = xtmp;
        xtmp = x[2];
        x[2] = x[3];
        x[3] = xtmp;
        xtmp = y[2];
        y[2] = y[3];
        y[3] = xtmp;
      }
      Vq = rtNaN;
      if ((!std::isnan(varargin_3)) && (!(varargin_3 > x[5])) &&
          (!(varargin_3 < x[0]))) {
        int high_i;
        int low_ip1;
        low_i = 1;
        low_ip1 = 2;
        high_i = 6;
        while (high_i > low_ip1) {
          int mid_i;
          mid_i = (low_i + high_i) >> 1;
          if (varargin_3 >= x[mid_i - 1]) {
            low_i = mid_i;
            low_ip1 = mid_i + 1;
          } else {
            high_i = mid_i;
          }
        }
        xtmp = x[low_i - 1];
        xtmp = (varargin_3 - xtmp) / (x[low_i] - xtmp);
        if (xtmp == 0.0) {
          Vq = y[low_i - 1];
        } else if (xtmp == 1.0) {
          Vq = y[low_i];
        } else {
          Vq = y[low_i - 1];
          if (!(Vq == y[low_i])) {
            Vq = (1.0 - xtmp) * Vq + xtmp * y[low_i];
          }
        }
      }
      exitg1 = 1;
    }
  } while (exitg1 == 0);
  return Vq;
}

//
// Arguments    : const double varargin_1[8]
//                const double varargin_2[8]
//                double varargin_3
// Return Type  : double
//
double c_interp1(const double varargin_1[8], const double varargin_2[8],
                 double varargin_3)
{
  double x[8];
  double y[8];
  double Vq;
  int low_i;
  std::copy(&varargin_2[0], &varargin_2[8], &y[0]);
  std::copy(&varargin_1[0], &varargin_1[8], &x[0]);
  low_i = 0;
  int exitg1;
  do {
    exitg1 = 0;
    if (low_i < 8) {
      if (std::isnan(varargin_1[low_i])) {
        exitg1 = 1;
      } else {
        low_i++;
      }
    } else {
      double xtmp;
      if (varargin_1[1] < varargin_1[0]) {
        xtmp = x[0];
        x[0] = x[7];
        x[7] = xtmp;
        xtmp = y[0];
        y[0] = y[7];
        y[7] = xtmp;
        xtmp = x[1];
        x[1] = x[6];
        x[6] = xtmp;
        xtmp = y[1];
        y[1] = y[6];
        y[6] = xtmp;
        xtmp = x[2];
        x[2] = x[5];
        x[5] = xtmp;
        xtmp = y[2];
        y[2] = y[5];
        y[5] = xtmp;
        xtmp = x[3];
        x[3] = x[4];
        x[4] = xtmp;
        xtmp = y[3];
        y[3] = y[4];
        y[4] = xtmp;
      }
      Vq = rtNaN;
      if ((!std::isnan(varargin_3)) && (!(varargin_3 > x[7])) &&
          (!(varargin_3 < x[0]))) {
        int high_i;
        int low_ip1;
        low_i = 1;
        low_ip1 = 2;
        high_i = 8;
        while (high_i > low_ip1) {
          int mid_i;
          mid_i = (low_i + high_i) >> 1;
          if (varargin_3 >= x[mid_i - 1]) {
            low_i = mid_i;
            low_ip1 = mid_i + 1;
          } else {
            high_i = mid_i;
          }
        }
        xtmp = x[low_i - 1];
        xtmp = (varargin_3 - xtmp) / (x[low_i] - xtmp);
        if (xtmp == 0.0) {
          Vq = y[low_i - 1];
        } else if (xtmp == 1.0) {
          Vq = y[low_i];
        } else {
          Vq = y[low_i - 1];
          if (!(Vq == y[low_i])) {
            Vq = (1.0 - xtmp) * Vq + xtmp * y[low_i];
          }
        }
      }
      exitg1 = 1;
    }
  } while (exitg1 == 0);
  return Vq;
}

//
// Arguments    : const double varargin_1[2]
//                const double varargin_2[2]
//                double varargin_3
// Return Type  : double
//
double d_interp1(const double varargin_1[2], const double varargin_2[2],
                 double varargin_3)
{
  double Vq;
  double r;
  double x_idx_1;
  double y_idx_0;
  double y_idx_1;
  int k;
  y_idx_0 = varargin_2[0];
  r = varargin_1[0];
  y_idx_1 = varargin_2[1];
  x_idx_1 = varargin_1[1];
  k = 0;
  int exitg1;
  do {
    exitg1 = 0;
    if (k < 2) {
      if (std::isnan(varargin_1[k])) {
        exitg1 = 1;
      } else {
        k++;
      }
    } else {
      if (varargin_1[1] < varargin_1[0]) {
        r = varargin_1[1];
        x_idx_1 = varargin_1[0];
        y_idx_0 = varargin_2[1];
        y_idx_1 = varargin_2[0];
      }
      Vq = rtNaN;
      if ((!std::isnan(varargin_3)) && (!(varargin_3 > x_idx_1)) &&
          (!(varargin_3 < r))) {
        r = (varargin_3 - r) / (x_idx_1 - r);
        if (r == 0.0) {
          Vq = y_idx_0;
        } else if (r == 1.0) {
          Vq = y_idx_1;
        } else if (y_idx_0 == y_idx_1) {
          Vq = y_idx_0;
        } else {
          Vq = (1.0 - r) * y_idx_0 + r * y_idx_1;
        }
      }
      exitg1 = 1;
    }
  } while (exitg1 == 0);
  return Vq;
}

//
// Arguments    : const double varargin_1[10]
//                const double varargin_2[10]
//                double varargin_3
// Return Type  : double
//
double interp1(const double varargin_1[10], const double varargin_2[10],
               double varargin_3)
{
  double x[10];
  double y[10];
  double Vq;
  int low_i;
  std::copy(&varargin_2[0], &varargin_2[10], &y[0]);
  std::copy(&varargin_1[0], &varargin_1[10], &x[0]);
  low_i = 0;
  int exitg1;
  do {
    exitg1 = 0;
    if (low_i < 10) {
      if (std::isnan(varargin_1[low_i])) {
        exitg1 = 1;
      } else {
        low_i++;
      }
    } else {
      double xtmp;
      if (varargin_1[1] < varargin_1[0]) {
        for (low_i = 0; low_i < 5; low_i++) {
          xtmp = x[low_i];
          x[low_i] = x[9 - low_i];
          x[9 - low_i] = xtmp;
          xtmp = y[low_i];
          y[low_i] = y[9 - low_i];
          y[9 - low_i] = xtmp;
        }
      }
      Vq = rtNaN;
      if ((!std::isnan(varargin_3)) && (!(varargin_3 > x[9])) &&
          (!(varargin_3 < x[0]))) {
        int high_i;
        int low_ip1;
        low_i = 1;
        low_ip1 = 2;
        high_i = 10;
        while (high_i > low_ip1) {
          int mid_i;
          mid_i = (low_i + high_i) >> 1;
          if (varargin_3 >= x[mid_i - 1]) {
            low_i = mid_i;
            low_ip1 = mid_i + 1;
          } else {
            high_i = mid_i;
          }
        }
        xtmp = x[low_i - 1];
        xtmp = (varargin_3 - xtmp) / (x[low_i] - xtmp);
        if (xtmp == 0.0) {
          Vq = y[low_i - 1];
        } else if (xtmp == 1.0) {
          Vq = y[low_i];
        } else {
          Vq = y[low_i - 1];
          if (!(Vq == y[low_i])) {
            Vq = (1.0 - xtmp) * Vq + xtmp * y[low_i];
          }
        }
      }
      exitg1 = 1;
    }
  } while (exitg1 == 0);
  return Vq;
}

} // namespace coder

//
// File trailer for interp1.cpp
//
// [EOF]
//
