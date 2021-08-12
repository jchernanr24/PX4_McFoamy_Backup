//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: _coder_McFoamy_FM_v3_info.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 05-Aug-2021 13:41:18
//

// Include Files
#include "_coder_McFoamy_FM_v3_info.h"
#include "emlrt.h"
#include "tmwtypes.h"

// Function Declarations
static const mxArray *emlrtMexFcnResolvedFunctionsInfo();

// Function Definitions
//
// Arguments    : void
// Return Type  : const mxArray *
//
static const mxArray *emlrtMexFcnResolvedFunctionsInfo()
{
  const mxArray *nameCaptureInfo;
  const char_T *data[7]{
      "789ced56cd6ed34010dea082b850222438219143cf5d4a0e28082105daa0a084364d4982"
      "a2caf8671d3bd9f5baf63a4db8c01bf00a7d0c8e3df2089c110f82ed"
      "c48e63b17214a781a61ec9194f3eef7e3363efb70b72d57a0e00b0ed5eafbe00f0f331f0"
      "eddec481fcd4df02f316c773537f3b1683f0ffadb9711efed2e5fb36",
      "8d656a30346293c010090a472a94e88668b093b18980856c8a8748f11155c7e84427a819"
      "0dde7b11a944a030f020effe8d86e441d321c0d2ec5986381af8fdf0"
      "ec3ba7dead05fbd1e2f4231fc3bb07a750a304c17ecffdd9a7b24390c16cd8d4898345a6"
      "5303b610a6b2cec6056459d42ac89a68f4103ca6cc870b0cd90cd6e5",
      "0a15c958a8d485617197cceaf89ab28e27097504b88a45d3444a59b754aa63a1563e1686"
      "7b42c52de1c8bd5098cfa794f9dce1e6334114ea4818cdeaff91924f"
      "e3f2cde32b7a8f096ddc25497dbcbf605d713f7bfeaeef7fd9873ec5baf802bb297c23ce"
      "7c8b7e970f397cf918fea1586cf4edc6bbb7bdb38e34aa575153156b",
      "913c8e127892f2009c785df36feafa3653d6b51d8be375053819ab8e2188d87c2a20554d"
      "5e0f49fc71e3f107b6aafdb69dc017e02b7a8ff36d7337dc75e9c8ef"
      "4c97af942fad2e3fe2f0e563382eb65b670dc71894dbe7c3cebec5ce4bac5fd91c5dbee0"
      "8c5fb48f2f38f3e76378b75aeb1c9cee10916151b228653b90518a25",
      "3a8288608875094e30484da61377d542f533b2a877468ee67b5df576d97371a6b77ff7b3"
      "e733bd5d075fa6b7ab997f53cfc169bf8f070975057828304d26623c"
      "d166cffed579f87249be60fe6e025f80af5a9fc3f64d37d84ca737836f5d3afd6ca0940e"
      "1b7b9262169fab3dad531a7cd4fbaf3747a72f38e3ffd773f175d7df",
      "65cfc79729f932fdbd5abec06e0a5fa6bfe9e6ff038f910220",
      ""};
  nameCaptureInfo = nullptr;
  emlrtNameCaptureMxArrayR2016a(&data[0], 7408U, &nameCaptureInfo);
  return nameCaptureInfo;
}

//
// Arguments    : void
// Return Type  : mxArray *
//
mxArray *emlrtMexFcnProperties()
{
  mxArray *xEntryPoints;
  mxArray *xInputs;
  mxArray *xResult;
  const char_T *epFieldName[6]{
      "Name",           "NumberOfInputs", "NumberOfOutputs",
      "ConstantInputs", "FullPath",       "TimeStamp"};
  const char_T *propFieldName[5]{"Version", "ResolvedFunctions", "EntryPoints",
                                 "CoverageInfo", "IsPolymorphic"};
  xEntryPoints =
      emlrtCreateStructMatrix(1, 1, 6, (const char_T **)&epFieldName[0]);
  xInputs = emlrtCreateLogicalMatrix(1, 10);
  emlrtSetField(xEntryPoints, 0, (const char_T *)"Name",
                emlrtMxCreateString((const char_T *)"McFoamy_FM_v3"));
  emlrtSetField(xEntryPoints, 0, (const char_T *)"NumberOfInputs",
                emlrtMxCreateDoubleScalar(10.0));
  emlrtSetField(xEntryPoints, 0, (const char_T *)"NumberOfOutputs",
                emlrtMxCreateDoubleScalar(3.0));
  emlrtSetField(xEntryPoints, 0, (const char_T *)"ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 0, (const char_T *)"FullPath",
      emlrtMxCreateString(
          (const char_T *)"/home/jgme/Documents/Simulation/Velocity error "
                          "change/Rotation test/McFoamy_FM_v3.m"));
  emlrtSetField(xEntryPoints, 0, (const char_T *)"TimeStamp",
                emlrtMxCreateDoubleScalar(738373.53810185182));
  xResult =
      emlrtCreateStructMatrix(1, 1, 5, (const char_T **)&propFieldName[0]);
  emlrtSetField(
      xResult, 0, (const char_T *)"Version",
      emlrtMxCreateString((const char_T *)"9.10.0.1710957 (R2021a) Update 4"));
  emlrtSetField(xResult, 0, (const char_T *)"ResolvedFunctions",
                (mxArray *)emlrtMexFcnResolvedFunctionsInfo());
  emlrtSetField(xResult, 0, (const char_T *)"EntryPoints", xEntryPoints);
  return xResult;
}

//
// File trailer for _coder_McFoamy_FM_v3_info.cpp
//
// [EOF]
//
