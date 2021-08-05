//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: _coder_McFoamy_FM_v3_api.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 05-Aug-2021 13:41:18
//

#ifndef _CODER_MCFOAMY_FM_V3_API_H
#define _CODER_MCFOAMY_FM_V3_API_H

// Include Files
#include "emlrt.h"
#include "tmwtypes.h"
#include <algorithm>
#include <cstring>

// Variable Declarations
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

// Function Declarations
void McFoamy_FM_v3(real_T Ail_def, real_T Elev_def, real_T Rud_def,
                   real_T Thr_com, real_T v_u, real_T v_v, real_T v_w,
                   real_T w_p, real_T w_q, real_T w_r, real_T F_total_body[3],
                   real_T M_total_body[3], real_T AER_lim[3]);

void McFoamy_FM_v3_api(const mxArray *const prhs[10], int32_T nlhs,
                       const mxArray *plhs[3]);

void McFoamy_FM_v3_atexit();

void McFoamy_FM_v3_initialize();

void McFoamy_FM_v3_terminate();

void McFoamy_FM_v3_xil_shutdown();

void McFoamy_FM_v3_xil_terminate();

#endif
//
// File trailer for _coder_McFoamy_FM_v3_api.h
//
// [EOF]
//
