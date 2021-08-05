//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: _coder_McFoamy_FM_v3_mex.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 05-Aug-2021 13:41:18
//

#ifndef _CODER_MCFOAMY_FM_V3_MEX_H
#define _CODER_MCFOAMY_FM_V3_MEX_H

// Include Files
#include "emlrt.h"
#include "mex.h"
#include "tmwtypes.h"

// Function Declarations
MEXFUNCTION_LINKAGE void mexFunction(int32_T nlhs, mxArray *plhs[],
                                     int32_T nrhs, const mxArray *prhs[]);

emlrtCTX mexFunctionCreateRootTLS();

void unsafe_McFoamy_FM_v3_mexFunction(int32_T nlhs, mxArray *plhs[3],
                                      int32_T nrhs, const mxArray *prhs[10]);

#endif
//
// File trailer for _coder_McFoamy_FM_v3_mex.h
//
// [EOF]
//
