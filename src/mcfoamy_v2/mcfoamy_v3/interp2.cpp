//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: interp2.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 05-Aug-2021 13:41:18
//

// Include Files
#include "interp2.h"
#include "rt_nonfinite.h"

// Function Definitions
//
// Arguments    : const double V[110]
//                double Xq
//                double Yq
//                double extrapval
// Return Type  : double
//
namespace coder {
double interp2_dispatch(const double V[110], double Xq, double Yq,
                        double extrapval)
{
  static const double b_dv[11]{
      0.0, 0.1, 0.2, 0.30000000000000004, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  double Vq;
  if ((Xq >= 0.0) && (Xq <= 1.0) && (Yq >= 0.0) && (Yq <= 1.5707963267948966)) {
    double d;
    double d1;
    double qx1;
    double rx;
    int b_low_i;
    int high_i;
    int low_i;
    int low_ip1;
    int mid_i;
    low_i = 1;
    low_ip1 = 2;
    high_i = 11;
    while (high_i > low_ip1) {
      mid_i = (low_i + high_i) >> 1;
      if (Xq >= b_dv[mid_i - 1]) {
        low_i = mid_i;
        low_ip1 = mid_i + 1;
      } else {
        high_i = mid_i;
      }
    }
    b_low_i = 1;
    low_ip1 = 2;
    high_i = 10;
    while (high_i > low_ip1) {
      mid_i = (b_low_i + high_i) >> 1;
      if (Yq >= 0.17453292519943295 * (static_cast<double>(mid_i) - 1.0)) {
        b_low_i = mid_i;
        low_ip1 = mid_i + 1;
      } else {
        high_i = mid_i;
      }
    }
    d = b_dv[low_i - 1];
    if (Xq == d) {
      low_ip1 = b_low_i + 10 * (low_i - 1);
      qx1 = V[low_ip1 - 1];
      Vq = V[low_ip1];
    } else if (Xq == b_dv[low_i]) {
      low_ip1 = b_low_i + 10 * low_i;
      qx1 = V[low_ip1 - 1];
      Vq = V[low_ip1];
    } else {
      rx = (Xq - d) / (b_dv[low_i] - d);
      low_ip1 = (b_low_i + 10 * (low_i - 1)) - 1;
      d = V[low_ip1];
      d1 = V[(b_low_i + 10 * low_i) - 1];
      if (d == d1) {
        qx1 = V[low_ip1];
      } else {
        qx1 = (1.0 - rx) * d + rx * d1;
      }
      low_ip1 = b_low_i + 10 * (low_i - 1);
      d = V[b_low_i + 10 * low_i];
      if (V[low_ip1] == d) {
        Vq = V[low_ip1];
      } else {
        Vq = (1.0 - rx) * V[low_ip1] + rx * d;
      }
    }
    d = 0.17453292519943295 * (static_cast<double>(b_low_i) - 1.0);
    if ((Yq == d) || (qx1 == Vq)) {
      Vq = qx1;
    } else {
      d1 = 0.17453292519943295 * static_cast<double>(b_low_i);
      if (!(Yq == d1)) {
        rx = (Yq - d) / (d1 - d);
        Vq = (1.0 - rx) * qx1 + rx * Vq;
      }
    }
  } else {
    Vq = extrapval;
  }
  return Vq;
}

} // namespace coder

//
// File trailer for interp2.cpp
//
// [EOF]
//
