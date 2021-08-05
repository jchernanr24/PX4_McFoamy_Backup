//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: flappedAirfoil_LAR_v1_FlatPlate.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 05-Aug-2021 13:41:18
//

#ifndef FLAPPEDAIRFOIL_LAR_V1_FLATPLATE_H
#define FLAPPEDAIRFOIL_LAR_V1_FLATPLATE_H

// Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
void b_flappedAirfoil_LAR_v1_FlatPlate(double alp, double cf, double c,
                                       double def, double Cd90, double *CN,
                                       double *CL, double *CD, double *CM);

void c_flappedAirfoil_LAR_v1_FlatPlate(double alp, double cf, double c,
                                       double def, double Cd90, double *CN,
                                       double *CL, double *CD, double *CM);

void flappedAirfoil_LAR_v1_FlatPlate(double alp, double cf, double c,
                                     double Cd90, double *CN, double *CL,
                                     double *CD, double *CM);

void flappedAirfoil_LAR_v1_FlatPlate(double alp, double cf, double c,
                                     double def, double Cd90, double *CN,
                                     double *CL, double *CD, double *CM);

#endif
//
// File trailer for flappedAirfoil_LAR_v1_FlatPlate.h
//
// [EOF]
//
