//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: sign.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 05-Aug-2021 13:41:18
//

// Include Files
#include "sign.h"
#include "rt_nonfinite.h"

// Function Definitions
//
// Arguments    : double x[4]
// Return Type  : void
//
namespace coder {
void b_sign(double x[4])
{
  double b_x;
  b_x = x[0];
  if (x[0] < 0.0) {
    b_x = -1.0;
  } else if (x[0] > 0.0) {
    b_x = 1.0;
  } else if (x[0] == 0.0) {
    b_x = 0.0;
  }
  x[0] = b_x;
  b_x = x[1];
  if (x[1] < 0.0) {
    b_x = -1.0;
  } else if (x[1] > 0.0) {
    b_x = 1.0;
  } else if (x[1] == 0.0) {
    b_x = 0.0;
  }
  x[1] = b_x;
  b_x = x[2];
  if (x[2] < 0.0) {
    b_x = -1.0;
  } else if (x[2] > 0.0) {
    b_x = 1.0;
  } else if (x[2] == 0.0) {
    b_x = 0.0;
  }
  x[2] = b_x;
  b_x = x[3];
  if (x[3] < 0.0) {
    b_x = -1.0;
  } else if (x[3] > 0.0) {
    b_x = 1.0;
  } else if (x[3] == 0.0) {
    b_x = 0.0;
  }
  x[3] = b_x;
}

} // namespace coder

//
// File trailer for sign.cpp
//
// [EOF]
//
