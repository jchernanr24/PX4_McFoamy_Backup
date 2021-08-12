//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: _coder_McFoamy_FM_v3_api.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 05-Aug-2021 13:41:18
//

// Include Files
#include "_coder_McFoamy_FM_v3_api.h"
#include "_coder_McFoamy_FM_v3_mex.h"

// Variable Definitions
emlrtCTX emlrtRootTLSGlobal{nullptr};

emlrtContext emlrtContextGlobal{
    true,                                                 // bFirstTime
    false,                                                // bInitialized
    131610U,                                              // fVersionInfo
    nullptr,                                              // fErrorFunction
    "McFoamy_FM_v3",                                      // fFunctionName
    nullptr,                                              // fRTCallStack
    false,                                                // bDebugMode
    {2045744189U, 2170104910U, 2743257031U, 4284093946U}, // fSigWrd
    nullptr                                               // fSigMem
};

// Function Declarations
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

static const mxArray *b_emlrt_marshallOut(const real_T u[3]);

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *Ail_def,
                               const char_T *identifier);

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId);

static const mxArray *emlrt_marshallOut(const real_T u[3]);

// Function Definitions
//
// Arguments    : const emlrtStack *sp
//                const mxArray *src
//                const emlrtMsgIdentifier *msgId
// Return Type  : real_T
//
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims{0};
  real_T ret;
  emlrtCheckBuiltInR2012b((emlrtCTX)sp, msgId, src, (const char_T *)"double",
                          false, 0U, (void *)&dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

//
// Arguments    : const real_T u[3]
// Return Type  : const mxArray *
//
static const mxArray *b_emlrt_marshallOut(const real_T u[3])
{
  static const int32_T i{0};
  static const int32_T i1{3};
  const mxArray *m;
  const mxArray *y;
  y = nullptr;
  m = emlrtCreateNumericArray(1, (const void *)&i, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m, &i1, 1);
  emlrtAssign(&y, m);
  return y;
}

//
// Arguments    : const emlrtStack *sp
//                const mxArray *Ail_def
//                const char_T *identifier
// Return Type  : real_T
//
static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *Ail_def,
                               const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = const_cast<const char_T *>(identifier);
  thisId.fParent = nullptr;
  thisId.bParentIsCell = false;
  y = emlrt_marshallIn(sp, emlrtAlias(Ail_def), &thisId);
  emlrtDestroyArray(&Ail_def);
  return y;
}

//
// Arguments    : const emlrtStack *sp
//                const mxArray *u
//                const emlrtMsgIdentifier *parentId
// Return Type  : real_T
//
static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = b_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

//
// Arguments    : const real_T u[3]
// Return Type  : const mxArray *
//
static const mxArray *emlrt_marshallOut(const real_T u[3])
{
  static const int32_T iv[2]{0, 0};
  static const int32_T iv1[2]{1, 3};
  const mxArray *m;
  const mxArray *y;
  y = nullptr;
  m = emlrtCreateNumericArray(2, (const void *)&iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m, &iv1[0], 2);
  emlrtAssign(&y, m);
  return y;
}

//
// Arguments    : const mxArray * const prhs[10]
//                int32_T nlhs
//                const mxArray *plhs[3]
// Return Type  : void
//
void McFoamy_FM_v3_api(const mxArray *const prhs[10], int32_T nlhs,
                       const mxArray *plhs[3])
{
  emlrtStack st{
      nullptr, // site
      nullptr, // tls
      nullptr  // prev
  };
  real_T(*AER_lim)[3];
  real_T(*F_total_body)[3];
  real_T(*M_total_body)[3];
  real_T Ail_def;
  real_T Elev_def;
  real_T Rud_def;
  real_T Thr_com;
  real_T v_u;
  real_T v_v;
  real_T v_w;
  real_T w_p;
  real_T w_q;
  real_T w_r;
  st.tls = emlrtRootTLSGlobal;
  F_total_body = (real_T(*)[3])mxMalloc(sizeof(real_T[3]));
  M_total_body = (real_T(*)[3])mxMalloc(sizeof(real_T[3]));
  AER_lim = (real_T(*)[3])mxMalloc(sizeof(real_T[3]));
  // Marshall function inputs
  Ail_def = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "Ail_def");
  Elev_def = emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "Elev_def");
  Rud_def = emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "Rud_def");
  Thr_com = emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "Thr_com");
  v_u = emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "v_u");
  v_v = emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "v_v");
  v_w = emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "v_w");
  w_p = emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "w_p");
  w_q = emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "w_q");
  w_r = emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "w_r");
  // Invoke the target function
  McFoamy_FM_v3(Ail_def, Elev_def, Rud_def, Thr_com, v_u, v_v, v_w, w_p, w_q,
                w_r, *F_total_body, *M_total_body, *AER_lim);
  // Marshall function outputs
  plhs[0] = emlrt_marshallOut(*F_total_body);
  if (nlhs > 1) {
    plhs[1] = emlrt_marshallOut(*M_total_body);
  }
  if (nlhs > 2) {
    plhs[2] = b_emlrt_marshallOut(*AER_lim);
  }
}

//
// Arguments    : void
// Return Type  : void
//
void McFoamy_FM_v3_atexit()
{
  emlrtStack st{
      nullptr, // site
      nullptr, // tls
      nullptr  // prev
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  McFoamy_FM_v3_xil_terminate();
  McFoamy_FM_v3_xil_shutdown();
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

//
// Arguments    : void
// Return Type  : void
//
void McFoamy_FM_v3_initialize()
{
  emlrtStack st{
      nullptr, // site
      nullptr, // tls
      nullptr  // prev
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, nullptr);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

//
// Arguments    : void
// Return Type  : void
//
void McFoamy_FM_v3_terminate()
{
  emlrtStack st{
      nullptr, // site
      nullptr, // tls
      nullptr  // prev
  };
  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

//
// File trailer for _coder_McFoamy_FM_v3_api.cpp
//
// [EOF]
//
