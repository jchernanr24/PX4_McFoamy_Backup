//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: flappedAirfoil_LAR_v1_FlatPlate.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 05-Aug-2021 13:41:18
//

// Include Files
#include "flappedAirfoil_LAR_v1_FlatPlate.h"
#include "fzero.h"
#include "interp1.h"
#include "rt_nonfinite.h"
#include <cmath>

// Variable Definitions
static const double dv[10]{0.1666, 0.333, 0.4, 0.5, 1.0,
                           1.25,   2.0,   3.0, 4.0, 6.0};

static const double dv1[10]{0.69813170079773179, 1.0471975511965976,
                            0.95993108859688125, 0.97738438111682457,
                            0.69813170079773179, 0.50614548307835561,
                            0.48869219055841229, 0.41887902047863906,
                            0.38397243543875248, 0.3490658503988659};

static const double dv2[8]{0.582, 0.536, 0.494, 0.456,
                           0.42,  0.388, 0.359, 0.333};

static const double dv3[6]{0.0, 0.2, 0.4, 0.6, 0.8, 1.0};

static const double dv4[6]{1.0, 0.83, 0.65, 0.3, 0.15, 0.0};

// Function Declarations
static double rt_powd_snf(double u0, double u1);

// Function Definitions
//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_powd_snf(double u0, double u1)
{
  double y;
  if (std::isnan(u0) || std::isnan(u1)) {
    y = rtNaN;
  } else {
    double d;
    double d1;
    d = std::abs(u0);
    d1 = std::abs(u1);
    if (std::isinf(u1)) {
      if (d == 1.0) {
        y = 1.0;
      } else if (d > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = rtNaN;
    } else {
      y = std::pow(u0, u1);
    }
  }
  return y;
}

//
// Arguments    : double alp
//                double cf
//                double c
//                double def
//                double Cd90
//                double *CN
//                double *CL
//                double *CD
//                double *CM
// Return Type  : void
//
void b_flappedAirfoil_LAR_v1_FlatPlate(double alp, double cf, double c,
                                       double def, double Cd90, double *CN,
                                       double *CL, double *CD, double *CM)
{
  double b_dv[8];
  double CN1[2];
  double b_LowAlpStart_P[2];
  double CT2;
  double HighAlpEnd_N;
  double HighAlpStart;
  double HighAlpStart_N;
  double LowAlpEnd_N;
  double LowAlpEnd_P;
  double LowAlpStart_N;
  double LowAlpStart_P;
  double b_HighAlpEnd_P;
  double b_gamma;
  double d;
  double delCLmax;
  double the_f;
  int HighAlpEnd_P;
  // ----------------- STALL ANGLES AND HIGH AOA ANGLES -----------------%
  //  AR_pre = [0.5,0.75,1,1.25,1.5,1.75,2, 3,4,6];
  //  alpStall_pre = [44,37,34,22,20,18,18, 18,17,17]*pi/180;
  //  alpLim_pre = [35,33,28,20,15,14,13, 13,13,13]*pi/180; % By Torres
  // alpLim_pre = [32,36,34,38,32,21,16,11,10,8]*pi/180; % Alp_Limit
  //  Alp_Limit
  HighAlpStart = coder::interp1(dv, dv1, 1.2420750265020692);
  // ----------------------- LIFT CURVE SLOPE ---------------------------%
  //  From McCormick
  //  -------- POTENTIAL LIFT AND VORTEX LIFT COEFFICIENTS --------------%
  // --------------------- FLAP DEFLECTION ------------------------------%
  d = cf / c;
  if (d == 1.0) {
    the_f = 0.0;
    LowAlpEnd_P = 0.37260508566838529;
    LowAlpEnd_N = -0.37260508566838529;
    LowAlpStart_P = 2.7689875679214078;
    LowAlpStart_N = -2.7689875679214078;
    b_HighAlpEnd_P = 3.1415926535897931 - HighAlpStart;
    HighAlpStart_N = -HighAlpStart;
    HighAlpEnd_N = -(3.1415926535897931 - HighAlpStart);
    alp += def;
    b_gamma = 0.0;
  } else {
    the_f = std::acos(2.0 * cf / c - 1.0);
    // eta_pre = [0.81*0.7 0.8*0.7 0.685*0.7 0.535 0.455 0.405 0.37 0.345];
    for (HighAlpEnd_P = 0; HighAlpEnd_P < 8; HighAlpEnd_P++) {
      b_dv[HighAlpEnd_P] =
          0.17453292519943295 * static_cast<double>(HighAlpEnd_P);
    }
    HighAlpEnd_N = std::abs(def);
    the_f = 1.7436202582132567 *
            (1.0 - (the_f - std::sin(the_f)) / 3.1415926535897931) *
            coder::c_interp1(b_dv, dv2, HighAlpEnd_N) * def;
    delCLmax = coder::b_interp1(dv3, dv4, d) * the_f;
    //  For Linear CL vs Alp curve
    //  alp0_eff1 = alp0 - delCL/CLAlp;
    //  CLmaxP1 = CLAlp*(alpStallP - alp0) + delCLmax;
    //  CLmaxN1 = CLAlp*(alpStallN - alp0) + delCLmax;
    //  alpStallP_eff1 = alp0_eff1 + CLmaxP1/CLAlp;
    //  alpStallN_eff1 = alp0_eff1 + CLmaxN1/CLAlp;
    //  For Non-Linear CL vs Alp curve
    // fun = @(x) Kp*sin(0-x)*(cos(0-x))^2 + Kv*abs(sin(0-x))*sin(0-x)*cos(0-x)
    // - delCL;
    the_f = coder::b_fzero(the_f);
    // fun = @(x) Kp*sin(x-alp0_eff)*(cos(x-alp0_eff))^2 +
    // Kv*abs(sin(x-alp0_eff))*sin(x-alp0_eff)*cos(x-alp0_eff) - CLmaxP;
    LowAlpEnd_P = coder::c_fzero(the_f, delCLmax + 0.93840864618177888);
    // alpStallP_eff = fzero(@myfun_alpStall_eff,x0,[],Kp,Kv,alp0_eff,CLmaxP);
    // fun = @(x) Kp*sin(x-alp0_eff)*(cos(x-alp0_eff))^2 +
    // Kv*abs(sin(x-alp0_eff))*sin(x-alp0_eff)*cos(x-alp0_eff) - CLmaxN;
    LowAlpEnd_N = coder::d_fzero(the_f, delCLmax + -0.93840864618177888);
    // alpStallN_eff = fzero(@myfun_alpStall_eff,x0,[],Kp,Kv,alp0_eff,CLmaxN);
    LowAlpStart_P = LowAlpEnd_N + 3.1415926535897931;
    LowAlpStart_N = LowAlpEnd_P - 3.1415926535897931;
    //  High AoA parameters using equivalent flat plate method
    CT2 = c - cf;
    b_gamma = std::asin(std::sin(def) * cf /
                        std::sqrt((CT2 * CT2 + cf * cf) +
                                  2.0 * CT2 * cf * std::cos(HighAlpEnd_N)));
    b_HighAlpEnd_P = 3.1415926535897931 - HighAlpStart;
    HighAlpStart_N = -HighAlpStart;
    HighAlpEnd_N = -(3.1415926535897931 - HighAlpStart);
  }
  //     %% Calculate Aerodynamic coefficients based on regime
  if ((alp >= -3.1415926535897931) && (alp <= LowAlpStart_N)) {
    double CL_tmp;
    double a_eff;
    double b_CL_tmp;
    a_eff = (alp - the_f) + 3.1415926535897931;
    CL_tmp = std::sin(a_eff);
    b_CL_tmp = std::abs(CL_tmp);
    the_f = std::cos(a_eff);
    *CL = 1.7436202582132567 * CL_tmp * (the_f * the_f) +
          3.1415926535897931 * b_CL_tmp * CL_tmp * the_f;
    delCLmax = std::sin(a_eff);
    *CD = std::abs(1.7436202582132567 * b_CL_tmp * CL_tmp * the_f +
                   3.1415926535897931 * rt_powd_snf(delCLmax, 3.0)) +
          0.02;
    *CN = -(1.7436202582132567 * delCLmax * std::abs(std::cos(a_eff)) +
            3.1415926535897931 * std::abs(delCLmax) * delCLmax);
    // CM = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
    *CM = -(0.25 - 0.175 * (1.0 - std::abs(a_eff - 3.1415926535897931) /
                                      1.5707963267948966)) *
          *CN;
    //  ------------------------------------------------- %
  } else if ((alp > LowAlpStart_N) && (alp < HighAlpEnd_N)) {
    double CL_tmp;
    double a;
    double a_eff;
    double b_CL_tmp;
    double y;
    //  Linear End
    a_eff = (LowAlpStart_N - the_f) + 3.1415926535897931;
    a = std::cos(a_eff);
    b_HighAlpEnd_P = std::sin(a_eff);
    LowAlpEnd_N = std::cos(a_eff);
    the_f = std::sin(a_eff);
    delCLmax = std::abs(the_f);
    HighAlpStart_N = -(1.7436202582132567 * the_f * LowAlpEnd_N +
                       3.1415926535897931 * delCLmax * the_f);
    // CM1 = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
    // CD1 = Cd0 + abs(CL1)*abs(tan(a_eff));
    // CN1 = -sqrt(CL1^2 + CD1^2);
    y = std::abs(a_eff - 3.1415926535897931);
    //  NonLinear End
    a_eff = HighAlpEnd_N + b_gamma;
    // - alp0_eff;
    //  Determining effective Cd90
    if ((d == 1.0) || (d == 0.0)) {
      b_gamma = Cd90;
    } else {
      if (HighAlpEnd_N < 0.0) {
        HighAlpEnd_P = -1;
      } else {
        HighAlpEnd_P = (HighAlpEnd_N > 0.0);
      }
      if (HighAlpEnd_P == 1) {
        b_CL_tmp = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (b_CL_tmp * b_CL_tmp) +
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      } else {
        b_CL_tmp = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (b_CL_tmp * b_CL_tmp) -
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      }
    }
    //  Hoerner Model with Lindenberg Correction
    b_CL_tmp = std::sin(std::abs(a_eff));
    b_gamma = -b_gamma *
              (1.0 / (0.44 * b_CL_tmp + 0.56) - 0.40999953366984909) * b_CL_tmp;
    b_CL_tmp = std::cos(a_eff);
    CT2 = 0.01 * b_CL_tmp;
    CL_tmp = std::sin(a_eff);
    //  Interpolation
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    CN1[0] = HighAlpStart_N;
    CN1[1] = b_gamma;
    *CN = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    CN1[0] = 1.7436202582132567 * b_HighAlpEnd_P * (a * a) +
             3.1415926535897931 * std::abs(b_HighAlpEnd_P) * b_HighAlpEnd_P *
                 LowAlpEnd_N;
    CN1[1] = b_gamma * b_CL_tmp - CT2 * CL_tmp;
    *CL = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    CN1[0] = std::abs(1.7436202582132567 * delCLmax * the_f * LowAlpEnd_N +
                      3.1415926535897931 * rt_powd_snf(the_f, 3.0)) +
             0.02;
    CN1[1] = b_gamma * CL_tmp + CT2 * b_CL_tmp;
    *CD = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    CN1[0] = -(0.25 - 0.175 * (1.0 - y / 1.5707963267948966)) * HighAlpStart_N;
    CN1[1] = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) *
             b_gamma;
    *CM = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    //  ------------------------------------------------- %
  } else if ((alp >= HighAlpEnd_N) && (alp <= HighAlpStart_N)) {
    double CL_tmp;
    double a_eff;
    //  Equivalent flat plate
    a_eff = alp + b_gamma;
    // - alp0_eff;
    //  Determining effective Cd90
    if ((d == 1.0) || (d == 0.0)) {
      b_gamma = Cd90;
    } else {
      if (alp < 0.0) {
        HighAlpEnd_P = -1;
      } else {
        HighAlpEnd_P = (alp > 0.0);
      }
      if (HighAlpEnd_P == 1) {
        double a;
        a = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (a * a) +
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      } else {
        CT2 = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (CT2 * CT2) - 0.00361111 * CT2) + Cd90;
      }
    }
    //  Hoerner Model
    the_f = std::sin(std::abs(a_eff));
    *CN =
        -b_gamma * (1.0 / (0.44 * the_f + 0.56) - 0.40999953366984909) * the_f;
    the_f = std::cos(a_eff);
    delCLmax = 0.01 * the_f;
    CL_tmp = std::sin(a_eff);
    *CL = *CN * the_f - delCLmax * CL_tmp;
    *CD = *CN * CL_tmp + delCLmax * the_f;
    *CM = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) * *CN;
    //  ------------------------------------------------- %
  } else if ((alp > HighAlpStart_N) && (alp < LowAlpEnd_N)) {
    double CL_tmp;
    double CN1_tmp;
    double a;
    double a_eff;
    double b_CL_tmp;
    double y;
    //  Linear End
    a_eff = LowAlpEnd_N - the_f;
    LowAlpStart_N = std::sin(a_eff);
    the_f = std::abs(LowAlpStart_N);
    delCLmax = std::cos(a_eff);
    a = std::cos(a_eff);
    HighAlpEnd_N = rt_powd_snf(std::sin(a_eff), 3.0);
    CL_tmp = std::sin(a_eff);
    y = std::abs(std::sin(a_eff));
    b_HighAlpEnd_P = std::sin(a_eff);
    //  NonLinear End
    a_eff = HighAlpStart_N + b_gamma;
    // - alp0_eff;
    if ((d == 1.0) || (d == 0.0)) {
      b_gamma = Cd90;
    } else {
      if (HighAlpStart_N < 0.0) {
        HighAlpEnd_P = -1;
      } else {
        HighAlpEnd_P = (HighAlpStart_N > 0.0);
      }
      if (HighAlpEnd_P == 1) {
        b_CL_tmp = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (b_CL_tmp * b_CL_tmp) +
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      } else {
        b_CL_tmp = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (b_CL_tmp * b_CL_tmp) -
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      }
    }
    b_CL_tmp = std::sin(std::abs(a_eff));
    b_gamma = -b_gamma *
              (1.0 / (0.44 * b_CL_tmp + 0.56) - 0.40999953366984909) * b_CL_tmp;
    b_CL_tmp = std::cos(a_eff);
    CT2 = 0.01 * b_CL_tmp;
    //  Interpolation
    b_LowAlpStart_P[0] = LowAlpEnd_N;
    b_LowAlpStart_P[1] = HighAlpStart_N;
    CN1[0] = 1.7436202582132567 * CL_tmp * delCLmax +
             3.1415926535897931 * y * b_HighAlpEnd_P;
    CN1[1] = b_gamma;
    *CN = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_N;
    b_LowAlpStart_P[1] = HighAlpStart_N;
    CN1[0] = 1.7436202582132567 * LowAlpStart_N * (a * a) +
             3.1415926535897931 * the_f * LowAlpStart_N * delCLmax;
    CN1_tmp = std::sin(a_eff);
    CN1[1] = b_gamma * b_CL_tmp - CT2 * CN1_tmp;
    *CL = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_N;
    b_LowAlpStart_P[1] = HighAlpStart_N;
    CN1[0] = std::abs(1.7436202582132567 * the_f * LowAlpStart_N * delCLmax +
                      3.1415926535897931 * HighAlpEnd_N) +
             0.02;
    CN1[1] = b_gamma * CN1_tmp + CT2 * std::cos(a_eff);
    *CD = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_N;
    b_LowAlpStart_P[1] = HighAlpStart_N;
    CN1[0] = -0.53407075111026481 * std::abs(b_HighAlpEnd_P) * b_HighAlpEnd_P;
    CN1[1] = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) *
             b_gamma;
    *CM = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    //  ------------------------------------------------- %
  } else if ((alp >= LowAlpEnd_N) && (alp <= LowAlpEnd_P)) {
    double a_eff;
    a_eff = alp - the_f;
    CT2 = std::cos(a_eff);
    *CL =
        1.7436202582132567 * std::sin(a_eff) * (CT2 * CT2) +
        3.1415926535897931 * std::abs(std::sin(a_eff)) * std::sin(a_eff) * CT2;
    delCLmax = std::sin(a_eff);
    the_f = std::abs(delCLmax);
    *CD = (std::abs(1.7436202582132567 * the_f * delCLmax * std::cos(a_eff)) +
           0.02) +
          3.1415926535897931 * std::abs(rt_powd_snf(delCLmax, 3.0));
    *CN = 1.7436202582132567 * delCLmax * std::cos(a_eff) +
          3.1415926535897931 * the_f * delCLmax;
    *CM = -0.53407075111026481 * std::abs(std::sin(a_eff)) * std::sin(a_eff);
    //  ------------------------------------------------- %
  } else if ((alp > LowAlpEnd_P) && (alp < HighAlpStart)) {
    double CL_tmp;
    double CN1_tmp;
    double a;
    double a_eff;
    double b_CL_tmp;
    double y;
    //  Linear End
    a_eff = LowAlpEnd_P - the_f;
    a = std::cos(a_eff);
    CL_tmp = std::sin(a_eff);
    y = std::abs(std::sin(a_eff));
    b_HighAlpEnd_P = std::sin(a_eff);
    LowAlpEnd_N = std::cos(a_eff);
    delCLmax = std::sin(a_eff);
    the_f = std::cos(a_eff);
    //  NonLinear End
    a_eff = HighAlpStart + b_gamma;
    // - alp0_eff;
    if ((d == 1.0) || (d == 0.0)) {
      b_gamma = Cd90;
    } else {
      if (HighAlpStart < 0.0) {
        HighAlpEnd_P = -1;
      } else {
        HighAlpEnd_P = (HighAlpStart > 0.0);
      }
      if (HighAlpEnd_P == 1) {
        b_CL_tmp = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (b_CL_tmp * b_CL_tmp) +
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      } else {
        b_CL_tmp = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (b_CL_tmp * b_CL_tmp) -
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      }
    }
    b_CL_tmp = std::sin(std::abs(a_eff));
    b_gamma = b_gamma * (1.0 / (0.44 * b_CL_tmp + 0.56) - 0.40999953366984909) *
              b_CL_tmp;
    b_CL_tmp = std::cos(a_eff);
    CT2 = 0.01 * b_CL_tmp;
    //  Interpolation
    b_LowAlpStart_P[0] = LowAlpEnd_P;
    b_LowAlpStart_P[1] = HighAlpStart;
    CN1_tmp = std::abs(delCLmax);
    CN1[0] = 1.7436202582132567 * delCLmax * the_f +
             3.1415926535897931 * CN1_tmp * delCLmax;
    CN1[1] = b_gamma;
    *CN = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_P;
    b_LowAlpStart_P[1] = HighAlpStart;
    CN1[0] = 1.7436202582132567 * CL_tmp * (a * a) +
             3.1415926535897931 * y * b_HighAlpEnd_P * LowAlpEnd_N;
    CN1[1] = b_gamma * b_CL_tmp - CT2 * std::sin(a_eff);
    *CL = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_P;
    b_LowAlpStart_P[1] = HighAlpStart;
    CN1[0] = std::abs(1.7436202582132567 * std::abs(b_HighAlpEnd_P) *
                          b_HighAlpEnd_P * LowAlpEnd_N +
                      3.1415926535897931 * rt_powd_snf(b_HighAlpEnd_P, 3.0)) +
             0.02;
    CN1[1] = b_gamma * std::sin(a_eff) + CT2 * b_CL_tmp;
    *CD = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_P;
    b_LowAlpStart_P[1] = HighAlpStart;
    CN1[0] = -0.53407075111026481 * CN1_tmp * delCLmax;
    CN1[1] = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) *
             b_gamma;
    *CM = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    //  ------------------------------------------------- %
  } else if ((alp >= HighAlpStart) && (alp <= b_HighAlpEnd_P)) {
    double CL_tmp;
    double a_eff;
    //  Equivalent flat plate
    a_eff = alp + b_gamma;
    // - alp0_eff;
    //  Determining effective Cd90
    if ((d == 1.0) || (d == 0.0)) {
      b_gamma = Cd90;
    } else {
      if (alp < 0.0) {
        HighAlpEnd_P = -1;
      } else {
        HighAlpEnd_P = (alp > 0.0);
      }
      if (HighAlpEnd_P == 1) {
        double a;
        a = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (a * a) +
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      } else {
        double a;
        a = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (a * a) -
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      }
    }
    //  Hoerner Model
    the_f = std::sin(std::abs(a_eff));
    *CN = b_gamma * (1.0 / (0.44 * the_f + 0.56) - 0.40999953366984909) * the_f;
    the_f = std::cos(a_eff);
    delCLmax = 0.01 * the_f;
    CL_tmp = std::sin(a_eff);
    *CL = *CN * the_f - delCLmax * CL_tmp;
    *CD = *CN * CL_tmp + delCLmax * the_f;
    *CM = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) * *CN;
    //  ------------------------------------------------- %
  } else if ((alp > b_HighAlpEnd_P) && (alp < LowAlpStart_P)) {
    double CL_tmp;
    double CN1_tmp;
    double a_eff;
    double b_CL_tmp;
    double y;
    //  Linear End
    a_eff = (LowAlpStart_P - the_f) - 3.1415926535897931;
    LowAlpStart_N = std::sin(a_eff);
    CT2 = std::cos(a_eff);
    y = std::abs(std::sin(a_eff));
    b_CL_tmp = std::sin(a_eff);
    LowAlpEnd_N = std::cos(a_eff);
    CN1_tmp = std::abs(b_CL_tmp);
    HighAlpStart_N = -(1.7436202582132567 * b_CL_tmp * LowAlpEnd_N +
                       3.1415926535897931 * CN1_tmp * b_CL_tmp);
    HighAlpEnd_N = std::abs(a_eff - 3.1415926535897931);
    // CM1 = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
    //  NonLinear End
    a_eff = b_HighAlpEnd_P + b_gamma;
    // - alp0_eff;
    //  Determining effective Cd90
    if ((d == 1.0) || (d == 0.0)) {
      b_gamma = Cd90;
    } else {
      if (b_HighAlpEnd_P < 0.0) {
        HighAlpEnd_P = -1;
      } else {
        HighAlpEnd_P = (b_HighAlpEnd_P > 0.0);
      }
      if (HighAlpEnd_P == 1) {
        double a;
        a = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (a * a) +
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      } else {
        double a;
        a = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (a * a) -
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      }
    }
    //  Hoerner Model
    b_gamma = b_gamma *
              (1.0 / (0.44 * std::sin(std::abs(a_eff)) + 0.56) -
               0.40999953366984909) *
              std::sin(std::abs(a_eff));
    delCLmax = 0.01 * std::cos(a_eff);
    CL_tmp = std::sin(a_eff);
    the_f = std::cos(a_eff);
    //  Interpolation
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = b_HighAlpEnd_P;
    CN1[0] = HighAlpStart_N;
    CN1[1] = b_gamma;
    *CN = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = b_HighAlpEnd_P;
    CN1[0] = 1.7436202582132567 * LowAlpStart_N * (CT2 * CT2) +
             3.1415926535897931 * y * LowAlpStart_N * CT2;
    CN1[1] = b_gamma * the_f - delCLmax * CL_tmp;
    *CL = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = b_HighAlpEnd_P;
    CN1[0] = std::abs(1.7436202582132567 * CN1_tmp * b_CL_tmp * LowAlpEnd_N +
                      3.1415926535897931 * rt_powd_snf(b_CL_tmp, 3.0)) +
             0.02;
    CN1[1] = b_gamma * CL_tmp + delCLmax * the_f;
    *CD = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = b_HighAlpEnd_P;
    CN1[0] = -(0.25 - 0.175 * (1.0 - HighAlpEnd_N / 1.5707963267948966)) *
             HighAlpStart_N;
    CN1[1] = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) *
             b_gamma;
    *CM = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    //  ------------------------------------------------- %
  } else {
    double CL_tmp;
    double a;
    double a_eff;
    double b_CL_tmp;
    a_eff = (alp - the_f) - 3.1415926535897931;
    a = std::cos(a_eff);
    CL_tmp = std::sin(a_eff);
    b_CL_tmp = std::cos(a_eff);
    the_f = std::abs(CL_tmp);
    delCLmax = 1.7436202582132567 * CL_tmp;
    HighAlpEnd_N = 3.1415926535897931 * the_f * CL_tmp;
    *CL = delCLmax * (a * a) + HighAlpEnd_N * b_CL_tmp;
    *CD = std::abs(1.7436202582132567 * the_f * CL_tmp * b_CL_tmp +
                   3.1415926535897931 * rt_powd_snf(CL_tmp, 3.0)) +
          0.02;
    *CN = -(delCLmax * b_CL_tmp + HighAlpEnd_N);
    *CM = -(0.25 - 0.175 * (1.0 - std::abs(a_eff - 3.1415926535897931) /
                                      1.5707963267948966)) *
          *CN;
    // CM = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
  }
}

//
// Arguments    : double alp
//                double cf
//                double c
//                double def
//                double Cd90
//                double *CN
//                double *CL
//                double *CD
//                double *CM
// Return Type  : void
//
void c_flappedAirfoil_LAR_v1_FlatPlate(double alp, double cf, double c,
                                       double def, double Cd90, double *CN,
                                       double *CL, double *CD, double *CM)
{
  double b_dv[8];
  double CN1[2];
  double b_LowAlpStart_P[2];
  double CT2;
  double HighAlpEnd_N;
  double HighAlpStart;
  double HighAlpStart_N;
  double LowAlpEnd_N;
  double LowAlpEnd_P;
  double LowAlpStart_N;
  double LowAlpStart_P;
  double b_HighAlpEnd_P;
  double b_gamma;
  double d;
  double delCLmax;
  double the_f;
  int HighAlpEnd_P;
  // ----------------- STALL ANGLES AND HIGH AOA ANGLES -----------------%
  //  AR_pre = [0.5,0.75,1,1.25,1.5,1.75,2, 3,4,6];
  //  alpStall_pre = [44,37,34,22,20,18,18, 18,17,17]*pi/180;
  //  alpLim_pre = [35,33,28,20,15,14,13, 13,13,13]*pi/180; % By Torres
  // alpLim_pre = [32,36,34,38,32,21,16,11,10,8]*pi/180; % Alp_Limit
  //  Alp_Limit
  HighAlpStart = coder::interp1(dv, dv1, 1.2984123165744528);
  // ----------------------- LIFT CURVE SLOPE ---------------------------%
  //  From McCormick
  //  -------- POTENTIAL LIFT AND VORTEX LIFT COEFFICIENTS --------------%
  // --------------------- FLAP DEFLECTION ------------------------------%
  d = cf / c;
  if (d == 1.0) {
    the_f = 0.0;
    LowAlpEnd_P = 0.36088611410052907;
    LowAlpEnd_N = -0.36088611410052907;
    LowAlpStart_P = 2.7807065394892643;
    LowAlpStart_N = -2.7807065394892643;
    b_HighAlpEnd_P = 3.1415926535897931 - HighAlpStart;
    HighAlpStart_N = -HighAlpStart;
    HighAlpEnd_N = -(3.1415926535897931 - HighAlpStart);
    alp += def;
    b_gamma = 0.0;
  } else {
    the_f = std::acos(2.0 * cf / c - 1.0);
    // eta_pre = [0.81*0.7 0.8*0.7 0.685*0.7 0.535 0.455 0.405 0.37 0.345];
    for (HighAlpEnd_P = 0; HighAlpEnd_P < 8; HighAlpEnd_P++) {
      b_dv[HighAlpEnd_P] =
          0.17453292519943295 * static_cast<double>(HighAlpEnd_P);
    }
    HighAlpEnd_N = std::abs(def);
    the_f = 1.8084579108343679 *
            (1.0 - (the_f - std::sin(the_f)) / 3.1415926535897931) *
            coder::c_interp1(b_dv, dv2, HighAlpEnd_N) * def;
    delCLmax = coder::b_interp1(dv3, dv4, d) * the_f;
    //  For Linear CL vs Alp curve
    //  alp0_eff1 = alp0 - delCL/CLAlp;
    //  CLmaxP1 = CLAlp*(alpStallP - alp0) + delCLmax;
    //  CLmaxN1 = CLAlp*(alpStallN - alp0) + delCLmax;
    //  alpStallP_eff1 = alp0_eff1 + CLmaxP1/CLAlp;
    //  alpStallN_eff1 = alp0_eff1 + CLmaxN1/CLAlp;
    //  For Non-Linear CL vs Alp curve
    // fun = @(x) Kp*sin(0-x)*(cos(0-x))^2 + Kv*abs(sin(0-x))*sin(0-x)*cos(0-x)
    // - delCL;
    the_f = coder::c_fzero(the_f);
    // fun = @(x) Kp*sin(x-alp0_eff)*(cos(x-alp0_eff))^2 +
    // Kv*abs(sin(x-alp0_eff))*sin(x-alp0_eff)*cos(x-alp0_eff) - CLmaxP;
    LowAlpEnd_P = coder::e_fzero(the_f, delCLmax + 0.92542259148461326);
    // alpStallP_eff = fzero(@myfun_alpStall_eff,x0,[],Kp,Kv,alp0_eff,CLmaxP);
    // fun = @(x) Kp*sin(x-alp0_eff)*(cos(x-alp0_eff))^2 +
    // Kv*abs(sin(x-alp0_eff))*sin(x-alp0_eff)*cos(x-alp0_eff) - CLmaxN;
    LowAlpEnd_N = coder::f_fzero(the_f, delCLmax + -0.92542259148461326);
    // alpStallN_eff = fzero(@myfun_alpStall_eff,x0,[],Kp,Kv,alp0_eff,CLmaxN);
    LowAlpStart_P = LowAlpEnd_N + 3.1415926535897931;
    LowAlpStart_N = LowAlpEnd_P - 3.1415926535897931;
    //  High AoA parameters using equivalent flat plate method
    CT2 = c - cf;
    b_gamma = std::asin(std::sin(def) * cf /
                        std::sqrt((CT2 * CT2 + cf * cf) +
                                  2.0 * CT2 * cf * std::cos(HighAlpEnd_N)));
    b_HighAlpEnd_P = 3.1415926535897931 - HighAlpStart;
    HighAlpStart_N = -HighAlpStart;
    HighAlpEnd_N = -(3.1415926535897931 - HighAlpStart);
  }
  //     %% Calculate Aerodynamic coefficients based on regime
  if ((alp >= -3.1415926535897931) && (alp <= LowAlpStart_N)) {
    double CL_tmp;
    double a_eff;
    double b_CL_tmp;
    a_eff = (alp - the_f) + 3.1415926535897931;
    CL_tmp = std::sin(a_eff);
    b_CL_tmp = std::abs(CL_tmp);
    the_f = std::cos(a_eff);
    *CL = 1.8084579108343679 * CL_tmp * (the_f * the_f) +
          3.1415926535897931 * b_CL_tmp * CL_tmp * the_f;
    delCLmax = std::sin(a_eff);
    *CD = std::abs(1.8084579108343679 * b_CL_tmp * CL_tmp * the_f +
                   3.1415926535897931 * rt_powd_snf(delCLmax, 3.0)) +
          0.02;
    *CN = -(1.8084579108343679 * delCLmax * std::abs(std::cos(a_eff)) +
            3.1415926535897931 * std::abs(delCLmax) * delCLmax);
    // CM = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
    *CM = -(0.25 - 0.175 * (1.0 - std::abs(a_eff - 3.1415926535897931) /
                                      1.5707963267948966)) *
          *CN;
    //  ------------------------------------------------- %
  } else if ((alp > LowAlpStart_N) && (alp < HighAlpEnd_N)) {
    double CL_tmp;
    double a;
    double a_eff;
    double b_CL_tmp;
    double y;
    //  Linear End
    a_eff = (LowAlpStart_N - the_f) + 3.1415926535897931;
    a = std::cos(a_eff);
    b_HighAlpEnd_P = std::sin(a_eff);
    LowAlpEnd_N = std::cos(a_eff);
    the_f = std::sin(a_eff);
    delCLmax = std::abs(the_f);
    HighAlpStart_N = -(1.8084579108343679 * the_f * LowAlpEnd_N +
                       3.1415926535897931 * delCLmax * the_f);
    // CM1 = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
    // CD1 = Cd0 + abs(CL1)*abs(tan(a_eff));
    // CN1 = -sqrt(CL1^2 + CD1^2);
    y = std::abs(a_eff - 3.1415926535897931);
    //  NonLinear End
    a_eff = HighAlpEnd_N + b_gamma;
    // - alp0_eff;
    //  Determining effective Cd90
    if ((d == 1.0) || (d == 0.0)) {
      b_gamma = Cd90;
    } else {
      if (HighAlpEnd_N < 0.0) {
        HighAlpEnd_P = -1;
      } else {
        HighAlpEnd_P = (HighAlpEnd_N > 0.0);
      }
      if (HighAlpEnd_P == 1) {
        b_CL_tmp = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (b_CL_tmp * b_CL_tmp) +
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      } else {
        b_CL_tmp = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (b_CL_tmp * b_CL_tmp) -
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      }
    }
    //  Hoerner Model with Lindenberg Correction
    b_CL_tmp = std::sin(std::abs(a_eff));
    b_gamma = -b_gamma *
              (1.0 / (0.44 * b_CL_tmp + 0.56) - 0.40999915549189647) * b_CL_tmp;
    b_CL_tmp = std::cos(a_eff);
    CT2 = 0.01 * b_CL_tmp;
    CL_tmp = std::sin(a_eff);
    //  Interpolation
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    CN1[0] = HighAlpStart_N;
    CN1[1] = b_gamma;
    *CN = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    CN1[0] = 1.8084579108343679 * b_HighAlpEnd_P * (a * a) +
             3.1415926535897931 * std::abs(b_HighAlpEnd_P) * b_HighAlpEnd_P *
                 LowAlpEnd_N;
    CN1[1] = b_gamma * b_CL_tmp - CT2 * CL_tmp;
    *CL = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    CN1[0] = std::abs(1.8084579108343679 * delCLmax * the_f * LowAlpEnd_N +
                      3.1415926535897931 * rt_powd_snf(the_f, 3.0)) +
             0.02;
    CN1[1] = b_gamma * CL_tmp + CT2 * b_CL_tmp;
    *CD = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    CN1[0] = -(0.25 - 0.175 * (1.0 - y / 1.5707963267948966)) * HighAlpStart_N;
    CN1[1] = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) *
             b_gamma;
    *CM = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    //  ------------------------------------------------- %
  } else if ((alp >= HighAlpEnd_N) && (alp <= HighAlpStart_N)) {
    double CL_tmp;
    double a_eff;
    //  Equivalent flat plate
    a_eff = alp + b_gamma;
    // - alp0_eff;
    //  Determining effective Cd90
    if ((d == 1.0) || (d == 0.0)) {
      b_gamma = Cd90;
    } else {
      if (alp < 0.0) {
        HighAlpEnd_P = -1;
      } else {
        HighAlpEnd_P = (alp > 0.0);
      }
      if (HighAlpEnd_P == 1) {
        double a;
        a = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (a * a) +
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      } else {
        CT2 = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (CT2 * CT2) - 0.00361111 * CT2) + Cd90;
      }
    }
    //  Hoerner Model
    the_f = std::sin(std::abs(a_eff));
    *CN =
        -b_gamma * (1.0 / (0.44 * the_f + 0.56) - 0.40999915549189647) * the_f;
    the_f = std::cos(a_eff);
    delCLmax = 0.01 * the_f;
    CL_tmp = std::sin(a_eff);
    *CL = *CN * the_f - delCLmax * CL_tmp;
    *CD = *CN * CL_tmp + delCLmax * the_f;
    *CM = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) * *CN;
    //  ------------------------------------------------- %
  } else if ((alp > HighAlpStart_N) && (alp < LowAlpEnd_N)) {
    double CL_tmp;
    double CN1_tmp;
    double a;
    double a_eff;
    double b_CL_tmp;
    double y;
    //  Linear End
    a_eff = LowAlpEnd_N - the_f;
    LowAlpStart_N = std::sin(a_eff);
    the_f = std::abs(LowAlpStart_N);
    delCLmax = std::cos(a_eff);
    a = std::cos(a_eff);
    HighAlpEnd_N = rt_powd_snf(std::sin(a_eff), 3.0);
    CL_tmp = std::sin(a_eff);
    y = std::abs(std::sin(a_eff));
    b_HighAlpEnd_P = std::sin(a_eff);
    //  NonLinear End
    a_eff = HighAlpStart_N + b_gamma;
    // - alp0_eff;
    if ((d == 1.0) || (d == 0.0)) {
      b_gamma = Cd90;
    } else {
      if (HighAlpStart_N < 0.0) {
        HighAlpEnd_P = -1;
      } else {
        HighAlpEnd_P = (HighAlpStart_N > 0.0);
      }
      if (HighAlpEnd_P == 1) {
        b_CL_tmp = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (b_CL_tmp * b_CL_tmp) +
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      } else {
        b_CL_tmp = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (b_CL_tmp * b_CL_tmp) -
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      }
    }
    b_CL_tmp = std::sin(std::abs(a_eff));
    b_gamma = -b_gamma *
              (1.0 / (0.44 * b_CL_tmp + 0.56) - 0.40999915549189647) * b_CL_tmp;
    b_CL_tmp = std::cos(a_eff);
    CT2 = 0.01 * b_CL_tmp;
    //  Interpolation
    b_LowAlpStart_P[0] = LowAlpEnd_N;
    b_LowAlpStart_P[1] = HighAlpStart_N;
    CN1[0] = 1.8084579108343679 * CL_tmp * delCLmax +
             3.1415926535897931 * y * b_HighAlpEnd_P;
    CN1[1] = b_gamma;
    *CN = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_N;
    b_LowAlpStart_P[1] = HighAlpStart_N;
    CN1[0] = 1.8084579108343679 * LowAlpStart_N * (a * a) +
             3.1415926535897931 * the_f * LowAlpStart_N * delCLmax;
    CN1_tmp = std::sin(a_eff);
    CN1[1] = b_gamma * b_CL_tmp - CT2 * CN1_tmp;
    *CL = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_N;
    b_LowAlpStart_P[1] = HighAlpStart_N;
    CN1[0] = std::abs(1.8084579108343679 * the_f * LowAlpStart_N * delCLmax +
                      3.1415926535897931 * HighAlpEnd_N) +
             0.02;
    CN1[1] = b_gamma * CN1_tmp + CT2 * std::cos(a_eff);
    *CD = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_N;
    b_LowAlpStart_P[1] = HighAlpStart_N;
    CN1[0] = -0.53407075111026481 * std::abs(b_HighAlpEnd_P) * b_HighAlpEnd_P;
    CN1[1] = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) *
             b_gamma;
    *CM = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    //  ------------------------------------------------- %
  } else if ((alp >= LowAlpEnd_N) && (alp <= LowAlpEnd_P)) {
    double a_eff;
    a_eff = alp - the_f;
    CT2 = std::cos(a_eff);
    *CL =
        1.8084579108343679 * std::sin(a_eff) * (CT2 * CT2) +
        3.1415926535897931 * std::abs(std::sin(a_eff)) * std::sin(a_eff) * CT2;
    delCLmax = std::sin(a_eff);
    the_f = std::abs(delCLmax);
    *CD = (std::abs(1.8084579108343679 * the_f * delCLmax * std::cos(a_eff)) +
           0.02) +
          3.1415926535897931 * std::abs(rt_powd_snf(delCLmax, 3.0));
    *CN = 1.8084579108343679 * delCLmax * std::cos(a_eff) +
          3.1415926535897931 * the_f * delCLmax;
    *CM = -0.53407075111026481 * std::abs(std::sin(a_eff)) * std::sin(a_eff);
    //  ------------------------------------------------- %
  } else if ((alp > LowAlpEnd_P) && (alp < HighAlpStart)) {
    double CL_tmp;
    double CN1_tmp;
    double a;
    double a_eff;
    double b_CL_tmp;
    double y;
    //  Linear End
    a_eff = LowAlpEnd_P - the_f;
    a = std::cos(a_eff);
    CL_tmp = std::sin(a_eff);
    y = std::abs(std::sin(a_eff));
    b_HighAlpEnd_P = std::sin(a_eff);
    LowAlpEnd_N = std::cos(a_eff);
    delCLmax = std::sin(a_eff);
    the_f = std::cos(a_eff);
    //  NonLinear End
    a_eff = HighAlpStart + b_gamma;
    // - alp0_eff;
    if ((d == 1.0) || (d == 0.0)) {
      b_gamma = Cd90;
    } else {
      if (HighAlpStart < 0.0) {
        HighAlpEnd_P = -1;
      } else {
        HighAlpEnd_P = (HighAlpStart > 0.0);
      }
      if (HighAlpEnd_P == 1) {
        b_CL_tmp = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (b_CL_tmp * b_CL_tmp) +
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      } else {
        b_CL_tmp = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (b_CL_tmp * b_CL_tmp) -
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      }
    }
    b_CL_tmp = std::sin(std::abs(a_eff));
    b_gamma = b_gamma * (1.0 / (0.44 * b_CL_tmp + 0.56) - 0.40999915549189647) *
              b_CL_tmp;
    b_CL_tmp = std::cos(a_eff);
    CT2 = 0.01 * b_CL_tmp;
    //  Interpolation
    b_LowAlpStart_P[0] = LowAlpEnd_P;
    b_LowAlpStart_P[1] = HighAlpStart;
    CN1_tmp = std::abs(delCLmax);
    CN1[0] = 1.8084579108343679 * delCLmax * the_f +
             3.1415926535897931 * CN1_tmp * delCLmax;
    CN1[1] = b_gamma;
    *CN = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_P;
    b_LowAlpStart_P[1] = HighAlpStart;
    CN1[0] = 1.8084579108343679 * CL_tmp * (a * a) +
             3.1415926535897931 * y * b_HighAlpEnd_P * LowAlpEnd_N;
    CN1[1] = b_gamma * b_CL_tmp - CT2 * std::sin(a_eff);
    *CL = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_P;
    b_LowAlpStart_P[1] = HighAlpStart;
    CN1[0] = std::abs(1.8084579108343679 * std::abs(b_HighAlpEnd_P) *
                          b_HighAlpEnd_P * LowAlpEnd_N +
                      3.1415926535897931 * rt_powd_snf(b_HighAlpEnd_P, 3.0)) +
             0.02;
    CN1[1] = b_gamma * std::sin(a_eff) + CT2 * b_CL_tmp;
    *CD = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_P;
    b_LowAlpStart_P[1] = HighAlpStart;
    CN1[0] = -0.53407075111026481 * CN1_tmp * delCLmax;
    CN1[1] = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) *
             b_gamma;
    *CM = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    //  ------------------------------------------------- %
  } else if ((alp >= HighAlpStart) && (alp <= b_HighAlpEnd_P)) {
    double CL_tmp;
    double a_eff;
    //  Equivalent flat plate
    a_eff = alp + b_gamma;
    // - alp0_eff;
    //  Determining effective Cd90
    if ((d == 1.0) || (d == 0.0)) {
      b_gamma = Cd90;
    } else {
      if (alp < 0.0) {
        HighAlpEnd_P = -1;
      } else {
        HighAlpEnd_P = (alp > 0.0);
      }
      if (HighAlpEnd_P == 1) {
        double a;
        a = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (a * a) +
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      } else {
        double a;
        a = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (a * a) -
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      }
    }
    //  Hoerner Model
    the_f = std::sin(std::abs(a_eff));
    *CN = b_gamma * (1.0 / (0.44 * the_f + 0.56) - 0.40999915549189647) * the_f;
    the_f = std::cos(a_eff);
    delCLmax = 0.01 * the_f;
    CL_tmp = std::sin(a_eff);
    *CL = *CN * the_f - delCLmax * CL_tmp;
    *CD = *CN * CL_tmp + delCLmax * the_f;
    *CM = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) * *CN;
    //  ------------------------------------------------- %
  } else if ((alp > b_HighAlpEnd_P) && (alp < LowAlpStart_P)) {
    double CL_tmp;
    double CN1_tmp;
    double a_eff;
    double b_CL_tmp;
    double y;
    //  Linear End
    a_eff = (LowAlpStart_P - the_f) - 3.1415926535897931;
    LowAlpStart_N = std::sin(a_eff);
    CT2 = std::cos(a_eff);
    y = std::abs(std::sin(a_eff));
    b_CL_tmp = std::sin(a_eff);
    LowAlpEnd_N = std::cos(a_eff);
    CN1_tmp = std::abs(b_CL_tmp);
    HighAlpStart_N = -(1.8084579108343679 * b_CL_tmp * LowAlpEnd_N +
                       3.1415926535897931 * CN1_tmp * b_CL_tmp);
    HighAlpEnd_N = std::abs(a_eff - 3.1415926535897931);
    // CM1 = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
    //  NonLinear End
    a_eff = b_HighAlpEnd_P + b_gamma;
    // - alp0_eff;
    //  Determining effective Cd90
    if ((d == 1.0) || (d == 0.0)) {
      b_gamma = Cd90;
    } else {
      if (b_HighAlpEnd_P < 0.0) {
        HighAlpEnd_P = -1;
      } else {
        HighAlpEnd_P = (b_HighAlpEnd_P > 0.0);
      }
      if (HighAlpEnd_P == 1) {
        double a;
        a = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (a * a) +
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      } else {
        double a;
        a = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (a * a) -
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      }
    }
    //  Hoerner Model
    b_gamma = b_gamma *
              (1.0 / (0.44 * std::sin(std::abs(a_eff)) + 0.56) -
               0.40999915549189647) *
              std::sin(std::abs(a_eff));
    delCLmax = 0.01 * std::cos(a_eff);
    CL_tmp = std::sin(a_eff);
    the_f = std::cos(a_eff);
    //  Interpolation
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = b_HighAlpEnd_P;
    CN1[0] = HighAlpStart_N;
    CN1[1] = b_gamma;
    *CN = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = b_HighAlpEnd_P;
    CN1[0] = 1.8084579108343679 * LowAlpStart_N * (CT2 * CT2) +
             3.1415926535897931 * y * LowAlpStart_N * CT2;
    CN1[1] = b_gamma * the_f - delCLmax * CL_tmp;
    *CL = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = b_HighAlpEnd_P;
    CN1[0] = std::abs(1.8084579108343679 * CN1_tmp * b_CL_tmp * LowAlpEnd_N +
                      3.1415926535897931 * rt_powd_snf(b_CL_tmp, 3.0)) +
             0.02;
    CN1[1] = b_gamma * CL_tmp + delCLmax * the_f;
    *CD = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = b_HighAlpEnd_P;
    CN1[0] = -(0.25 - 0.175 * (1.0 - HighAlpEnd_N / 1.5707963267948966)) *
             HighAlpStart_N;
    CN1[1] = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) *
             b_gamma;
    *CM = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    //  ------------------------------------------------- %
  } else {
    double CL_tmp;
    double a;
    double a_eff;
    double b_CL_tmp;
    a_eff = (alp - the_f) - 3.1415926535897931;
    a = std::cos(a_eff);
    CL_tmp = std::sin(a_eff);
    b_CL_tmp = std::cos(a_eff);
    the_f = std::abs(CL_tmp);
    delCLmax = 1.8084579108343679 * CL_tmp;
    HighAlpEnd_N = 3.1415926535897931 * the_f * CL_tmp;
    *CL = delCLmax * (a * a) + HighAlpEnd_N * b_CL_tmp;
    *CD = std::abs(1.8084579108343679 * the_f * CL_tmp * b_CL_tmp +
                   3.1415926535897931 * rt_powd_snf(CL_tmp, 3.0)) +
          0.02;
    *CN = -(delCLmax * b_CL_tmp + HighAlpEnd_N);
    *CM = -(0.25 - 0.175 * (1.0 - std::abs(a_eff - 3.1415926535897931) /
                                      1.5707963267948966)) *
          *CN;
    // CM = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
  }
}

//
// Arguments    : double alp
//                double cf
//                double c
//                double Cd90
//                double *CN
//                double *CL
//                double *CD
//                double *CM
// Return Type  : void
//
void flappedAirfoil_LAR_v1_FlatPlate(double alp, double cf, double c,
                                     double Cd90, double *CN, double *CL,
                                     double *CD, double *CM)
{
  double b_dv[8];
  double CN1[2];
  double b_LowAlpStart_P[2];
  double HighAlpEnd_N;
  double HighAlpEnd_P;
  double HighAlpStart;
  double HighAlpStart_N;
  double LowAlpEnd_N;
  double LowAlpEnd_P;
  double LowAlpStart_N;
  double LowAlpStart_P;
  double a;
  double a_tmp;
  double b_gamma;
  double delCLmax;
  double the_f;
  // ----------------- STALL ANGLES AND HIGH AOA ANGLES -----------------%
  //  AR_pre = [0.5,0.75,1,1.25,1.5,1.75,2, 3,4,6];
  //  alpStall_pre = [44,37,34,22,20,18,18, 18,17,17]*pi/180;
  //  alpLim_pre = [35,33,28,20,15,14,13, 13,13,13]*pi/180; % By Torres
  // alpLim_pre = [32,36,34,38,32,21,16,11,10,8]*pi/180; % Alp_Limit
  //  Alp_Limit
  HighAlpStart = coder::interp1(dv, dv1, 0.22793785354755264);
  // ----------------------- LIFT CURVE SLOPE ---------------------------%
  //  From McCormick
  //  -------- POTENTIAL LIFT AND VORTEX LIFT COEFFICIENTS --------------%
  // --------------------- FLAP DEFLECTION ------------------------------%
  a = cf / c;
  if (a == 1.0) {
    the_f = 0.0;
    LowAlpEnd_P = 0.58423967555431222;
    LowAlpEnd_N = -0.58423967555431222;
    LowAlpStart_P = 2.5573529780354809;
    LowAlpStart_N = -2.5573529780354809;
    HighAlpEnd_P = 3.1415926535897931 - HighAlpStart;
    HighAlpStart_N = -HighAlpStart;
    HighAlpEnd_N = -(3.1415926535897931 - HighAlpStart);
    b_gamma = 0.0;
  } else {
    the_f = std::acos(2.0 * cf / c - 1.0);
    // eta_pre = [0.81*0.7 0.8*0.7 0.685*0.7 0.535 0.455 0.405 0.37 0.345];
    for (int i{0}; i < 8; i++) {
      b_dv[i] = 0.17453292519943295 * static_cast<double>(i);
    }
    the_f = 0.35596863975551224 *
            (1.0 - (the_f - std::sin(the_f)) / 3.1415926535897931) *
            coder::c_interp1(b_dv, dv2, 0.0) * 0.0;
    delCLmax = coder::b_interp1(dv3, dv4, a) * the_f;
    //  For Linear CL vs Alp curve
    //  alp0_eff1 = alp0 - delCL/CLAlp;
    //  CLmaxP1 = CLAlp*(alpStallP - alp0) + delCLmax;
    //  CLmaxN1 = CLAlp*(alpStallN - alp0) + delCLmax;
    //  alpStallP_eff1 = alp0_eff1 + CLmaxP1/CLAlp;
    //  alpStallN_eff1 = alp0_eff1 + CLmaxN1/CLAlp;
    //  For Non-Linear CL vs Alp curve
    // fun = @(x) Kp*sin(0-x)*(cos(0-x))^2 + Kv*abs(sin(0-x))*sin(0-x)*cos(0-x)
    // - delCL;
    the_f = coder::d_fzero(the_f);
    // fun = @(x) Kp*sin(x-alp0_eff)*(cos(x-alp0_eff))^2 +
    // Kv*abs(sin(x-alp0_eff))*sin(x-alp0_eff)*cos(x-alp0_eff) - CLmaxP;
    LowAlpEnd_P = coder::g_fzero(the_f, delCLmax + 0.93382899733155733);
    // alpStallP_eff = fzero(@myfun_alpStall_eff,x0,[],Kp,Kv,alp0_eff,CLmaxP);
    // fun = @(x) Kp*sin(x-alp0_eff)*(cos(x-alp0_eff))^2 +
    // Kv*abs(sin(x-alp0_eff))*sin(x-alp0_eff)*cos(x-alp0_eff) - CLmaxN;
    LowAlpEnd_N = coder::h_fzero(the_f, delCLmax + -0.93382899733155733);
    // alpStallN_eff = fzero(@myfun_alpStall_eff,x0,[],Kp,Kv,alp0_eff,CLmaxN);
    LowAlpStart_P = LowAlpEnd_N + 3.1415926535897931;
    LowAlpStart_N = LowAlpEnd_P - 3.1415926535897931;
    //  High AoA parameters using equivalent flat plate method
    a_tmp = c - cf;
    b_gamma = std::asin(
        0.0 * cf / std::sqrt((a_tmp * a_tmp + cf * cf) + 2.0 * a_tmp * cf));
    HighAlpEnd_P = 3.1415926535897931 - HighAlpStart;
    HighAlpStart_N = -HighAlpStart;
    HighAlpEnd_N = -(3.1415926535897931 - HighAlpStart);
  }
  //     %% Calculate Aerodynamic coefficients based on regime
  if ((alp >= -3.1415926535897931) && (alp <= LowAlpStart_N)) {
    double a_eff;
    a_eff = (alp - the_f) + 3.1415926535897931;
    a_tmp = std::sin(a_eff);
    HighAlpEnd_N = std::abs(a_tmp);
    the_f = std::cos(a_eff);
    *CL = 0.35596863975551224 * a_tmp * (the_f * the_f) +
          3.1415926535897931 * HighAlpEnd_N * a_tmp * the_f;
    delCLmax = std::sin(a_eff);
    *CD = std::abs(0.35596863975551224 * HighAlpEnd_N * a_tmp * the_f +
                   3.1415926535897931 * rt_powd_snf(delCLmax, 3.0)) +
          0.02;
    *CN = -(0.35596863975551224 * delCLmax * std::abs(std::cos(a_eff)) +
            3.1415926535897931 * std::abs(delCLmax) * delCLmax);
    // CM = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
    *CM = -(0.25 - 0.175 * (1.0 - std::abs(a_eff - 3.1415926535897931) /
                                      1.5707963267948966)) *
          *CN;
    //  ------------------------------------------------- %
  } else if ((alp > LowAlpStart_N) && (alp < HighAlpEnd_N)) {
    double a_eff;
    double y;
    //  Linear End
    a_eff = (LowAlpStart_N - the_f) + 3.1415926535897931;
    a = std::cos(a_eff);
    HighAlpEnd_P = std::sin(a_eff);
    LowAlpEnd_N = std::cos(a_eff);
    the_f = std::sin(a_eff);
    delCLmax = std::abs(the_f);
    LowAlpEnd_P = -(0.35596863975551224 * the_f * LowAlpEnd_N +
                    3.1415926535897931 * delCLmax * the_f);
    // CM1 = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
    // CD1 = Cd0 + abs(CL1)*abs(tan(a_eff));
    // CN1 = -sqrt(CL1^2 + CD1^2);
    y = std::abs(a_eff - 3.1415926535897931);
    //  NonLinear End
    a_eff = HighAlpEnd_N + b_gamma;
    // - alp0_eff;
    //  Determining effective Cd90
    //  Hoerner Model with Lindenberg Correction
    a_tmp = std::sin(std::abs(a_eff));
    b_gamma = -Cd90 * (1.0 / (0.44 * a_tmp + 0.56) - 0.41) * a_tmp;
    a_tmp = std::cos(a_eff);
    LowAlpStart_P = 0.01 * a_tmp;
    HighAlpStart_N = std::sin(a_eff);
    //  Interpolation
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    CN1[0] = LowAlpEnd_P;
    CN1[1] = b_gamma;
    *CN = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    CN1[0] = 0.35596863975551224 * HighAlpEnd_P * (a * a) +
             3.1415926535897931 * std::abs(HighAlpEnd_P) * HighAlpEnd_P *
                 LowAlpEnd_N;
    CN1[1] = b_gamma * a_tmp - LowAlpStart_P * HighAlpStart_N;
    *CL = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    CN1[0] = std::abs(0.35596863975551224 * delCLmax * the_f * LowAlpEnd_N +
                      3.1415926535897931 * rt_powd_snf(the_f, 3.0)) +
             0.02;
    CN1[1] = b_gamma * HighAlpStart_N + LowAlpStart_P * a_tmp;
    *CD = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    CN1[0] = -(0.25 - 0.175 * (1.0 - y / 1.5707963267948966)) * LowAlpEnd_P;
    CN1[1] = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) *
             b_gamma;
    *CM = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    //  ------------------------------------------------- %
  } else if ((alp >= HighAlpEnd_N) && (alp <= HighAlpStart_N)) {
    double a_eff;
    //  Equivalent flat plate
    a_eff = alp + b_gamma;
    // - alp0_eff;
    //  Determining effective Cd90
    //  Hoerner Model
    the_f = std::sin(std::abs(a_eff));
    *CN = -Cd90 * (1.0 / (0.44 * the_f + 0.56) - 0.41) * the_f;
    the_f = std::cos(a_eff);
    LowAlpStart_N = 0.01 * the_f;
    a_tmp = std::sin(a_eff);
    *CL = *CN * the_f - LowAlpStart_N * a_tmp;
    *CD = *CN * a_tmp + LowAlpStart_N * the_f;
    *CM = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) * *CN;
    //  ------------------------------------------------- %
  } else if ((alp > HighAlpStart_N) && (alp < LowAlpEnd_N)) {
    double CN1_tmp;
    double a_eff;
    double y;
    //  Linear End
    a_eff = LowAlpEnd_N - the_f;
    HighAlpStart = std::sin(a_eff);
    the_f = std::abs(HighAlpStart);
    delCLmax = std::cos(a_eff);
    a = std::cos(a_eff);
    HighAlpEnd_N = rt_powd_snf(std::sin(a_eff), 3.0);
    LowAlpStart_N = std::sin(a_eff);
    y = std::abs(std::sin(a_eff));
    HighAlpEnd_P = std::sin(a_eff);
    //  NonLinear End
    a_eff = HighAlpStart_N + b_gamma;
    // - alp0_eff;
    a_tmp = std::sin(std::abs(a_eff));
    b_gamma = -Cd90 * (1.0 / (0.44 * a_tmp + 0.56) - 0.41) * a_tmp;
    a_tmp = std::cos(a_eff);
    LowAlpStart_P = 0.01 * a_tmp;
    //  Interpolation
    b_LowAlpStart_P[0] = LowAlpEnd_N;
    b_LowAlpStart_P[1] = HighAlpStart_N;
    CN1[0] = 0.35596863975551224 * LowAlpStart_N * delCLmax +
             3.1415926535897931 * y * HighAlpEnd_P;
    CN1[1] = b_gamma;
    *CN = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_N;
    b_LowAlpStart_P[1] = HighAlpStart_N;
    CN1[0] = 0.35596863975551224 * HighAlpStart * (a * a) +
             3.1415926535897931 * the_f * HighAlpStart * delCLmax;
    CN1_tmp = std::sin(a_eff);
    CN1[1] = b_gamma * a_tmp - LowAlpStart_P * CN1_tmp;
    *CL = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_N;
    b_LowAlpStart_P[1] = HighAlpStart_N;
    CN1[0] = std::abs(0.35596863975551224 * the_f * HighAlpStart * delCLmax +
                      3.1415926535897931 * HighAlpEnd_N) +
             0.02;
    CN1[1] = b_gamma * CN1_tmp + LowAlpStart_P * std::cos(a_eff);
    *CD = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_N;
    b_LowAlpStart_P[1] = HighAlpStart_N;
    CN1[0] = -0.53407075111026481 * std::abs(HighAlpEnd_P) * HighAlpEnd_P;
    CN1[1] = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) *
             b_gamma;
    *CM = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    //  ------------------------------------------------- %
  } else if ((alp >= LowAlpEnd_N) && (alp <= LowAlpEnd_P)) {
    double a_eff;
    a_eff = alp - the_f;
    a_tmp = std::cos(a_eff);
    *CL = 0.35596863975551224 * std::sin(a_eff) * (a_tmp * a_tmp) +
          3.1415926535897931 * std::abs(std::sin(a_eff)) * std::sin(a_eff) *
              a_tmp;
    delCLmax = std::sin(a_eff);
    the_f = std::abs(delCLmax);
    *CD = (std::abs(0.35596863975551224 * the_f * delCLmax * std::cos(a_eff)) +
           0.02) +
          3.1415926535897931 * std::abs(rt_powd_snf(delCLmax, 3.0));
    *CN = 0.35596863975551224 * delCLmax * std::cos(a_eff) +
          3.1415926535897931 * the_f * delCLmax;
    *CM = -0.53407075111026481 * std::abs(std::sin(a_eff)) * std::sin(a_eff);
    //  ------------------------------------------------- %
  } else if ((alp > LowAlpEnd_P) && (alp < HighAlpStart)) {
    double CN1_tmp;
    double a_eff;
    double y;
    //  Linear End
    a_eff = LowAlpEnd_P - the_f;
    a = std::cos(a_eff);
    LowAlpStart_N = std::sin(a_eff);
    y = std::abs(std::sin(a_eff));
    HighAlpEnd_P = std::sin(a_eff);
    LowAlpEnd_N = std::cos(a_eff);
    delCLmax = std::sin(a_eff);
    the_f = std::cos(a_eff);
    //  NonLinear End
    a_eff = HighAlpStart + b_gamma;
    // - alp0_eff;
    a_tmp = std::sin(std::abs(a_eff));
    b_gamma = Cd90 * (1.0 / (0.44 * a_tmp + 0.56) - 0.41) * a_tmp;
    a_tmp = std::cos(a_eff);
    LowAlpStart_P = 0.01 * a_tmp;
    //  Interpolation
    b_LowAlpStart_P[0] = LowAlpEnd_P;
    b_LowAlpStart_P[1] = HighAlpStart;
    CN1_tmp = std::abs(delCLmax);
    CN1[0] = 0.35596863975551224 * delCLmax * the_f +
             3.1415926535897931 * CN1_tmp * delCLmax;
    CN1[1] = b_gamma;
    *CN = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_P;
    b_LowAlpStart_P[1] = HighAlpStart;
    CN1[0] = 0.35596863975551224 * LowAlpStart_N * (a * a) +
             3.1415926535897931 * y * HighAlpEnd_P * LowAlpEnd_N;
    CN1[1] = b_gamma * a_tmp - LowAlpStart_P * std::sin(a_eff);
    *CL = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_P;
    b_LowAlpStart_P[1] = HighAlpStart;
    CN1[0] = std::abs(0.35596863975551224 * std::abs(HighAlpEnd_P) *
                          HighAlpEnd_P * LowAlpEnd_N +
                      3.1415926535897931 * rt_powd_snf(HighAlpEnd_P, 3.0)) +
             0.02;
    CN1[1] = b_gamma * std::sin(a_eff) + LowAlpStart_P * a_tmp;
    *CD = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_P;
    b_LowAlpStart_P[1] = HighAlpStart;
    CN1[0] = -0.53407075111026481 * CN1_tmp * delCLmax;
    CN1[1] = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) *
             b_gamma;
    *CM = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    //  ------------------------------------------------- %
  } else if ((alp >= HighAlpStart) && (alp <= HighAlpEnd_P)) {
    double a_eff;
    //  Equivalent flat plate
    a_eff = alp + b_gamma;
    // - alp0_eff;
    //  Determining effective Cd90
    //  Hoerner Model
    the_f = std::sin(std::abs(a_eff));
    *CN = Cd90 * (1.0 / (0.44 * the_f + 0.56) - 0.41) * the_f;
    the_f = std::cos(a_eff);
    LowAlpStart_N = 0.01 * the_f;
    a_tmp = std::sin(a_eff);
    *CL = *CN * the_f - LowAlpStart_N * a_tmp;
    *CD = *CN * a_tmp + LowAlpStart_N * the_f;
    *CM = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) * *CN;
    //  ------------------------------------------------- %
  } else if ((alp > HighAlpEnd_P) && (alp < LowAlpStart_P)) {
    double CN1_tmp;
    double a_eff;
    double y;
    //  Linear End
    a_eff = (LowAlpStart_P - the_f) - 3.1415926535897931;
    HighAlpStart = std::sin(a_eff);
    a_tmp = std::cos(a_eff);
    y = std::abs(std::sin(a_eff));
    delCLmax = std::sin(a_eff);
    LowAlpEnd_N = std::cos(a_eff);
    CN1_tmp = std::abs(delCLmax);
    LowAlpEnd_P = -(0.35596863975551224 * delCLmax * LowAlpEnd_N +
                    3.1415926535897931 * CN1_tmp * delCLmax);
    HighAlpEnd_N = std::abs(a_eff - 3.1415926535897931);
    // CM1 = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
    //  NonLinear End
    a_eff = HighAlpEnd_P + b_gamma;
    // - alp0_eff;
    //  Determining effective Cd90
    //  Hoerner Model
    b_gamma = Cd90 * (1.0 / (0.44 * std::sin(std::abs(a_eff)) + 0.56) - 0.41) *
              std::sin(std::abs(a_eff));
    LowAlpStart_N = 0.01 * std::cos(a_eff);
    HighAlpStart_N = std::sin(a_eff);
    the_f = std::cos(a_eff);
    //  Interpolation
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = HighAlpEnd_P;
    CN1[0] = LowAlpEnd_P;
    CN1[1] = b_gamma;
    *CN = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = HighAlpEnd_P;
    CN1[0] = 0.35596863975551224 * HighAlpStart * (a_tmp * a_tmp) +
             3.1415926535897931 * y * HighAlpStart * a_tmp;
    CN1[1] = b_gamma * the_f - LowAlpStart_N * HighAlpStart_N;
    *CL = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = HighAlpEnd_P;
    CN1[0] = std::abs(0.35596863975551224 * CN1_tmp * delCLmax * LowAlpEnd_N +
                      3.1415926535897931 * rt_powd_snf(delCLmax, 3.0)) +
             0.02;
    CN1[1] = b_gamma * HighAlpStart_N + LowAlpStart_N * the_f;
    *CD = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = HighAlpEnd_P;
    CN1[0] = -(0.25 - 0.175 * (1.0 - HighAlpEnd_N / 1.5707963267948966)) *
             LowAlpEnd_P;
    CN1[1] = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) *
             b_gamma;
    *CM = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    //  ------------------------------------------------- %
  } else {
    double a_eff;
    a_eff = (alp - the_f) - 3.1415926535897931;
    a = std::cos(a_eff);
    a_tmp = std::sin(a_eff);
    HighAlpEnd_N = std::cos(a_eff);
    the_f = std::abs(a_tmp);
    delCLmax = 0.35596863975551224 * a_tmp;
    LowAlpStart_N = 3.1415926535897931 * the_f * a_tmp;
    *CL = delCLmax * (a * a) + LowAlpStart_N * HighAlpEnd_N;
    *CD = std::abs(0.35596863975551224 * the_f * a_tmp * HighAlpEnd_N +
                   3.1415926535897931 * rt_powd_snf(a_tmp, 3.0)) +
          0.02;
    *CN = -(delCLmax * HighAlpEnd_N + LowAlpStart_N);
    *CM = -(0.25 - 0.175 * (1.0 - std::abs(a_eff - 3.1415926535897931) /
                                      1.5707963267948966)) *
          *CN;
    // CM = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
  }
}

//
// Arguments    : double alp
//                double cf
//                double c
//                double def
//                double Cd90
//                double *CN
//                double *CL
//                double *CD
//                double *CM
// Return Type  : void
//
void flappedAirfoil_LAR_v1_FlatPlate(double alp, double cf, double c,
                                     double def, double Cd90, double *CN,
                                     double *CL, double *CD, double *CM)
{
  double b_dv[8];
  double CN1[2];
  double b_LowAlpStart_P[2];
  double CT2;
  double HighAlpEnd_N;
  double HighAlpStart;
  double HighAlpStart_N;
  double LowAlpEnd_N;
  double LowAlpEnd_P;
  double LowAlpStart_N;
  double LowAlpStart_P;
  double b_HighAlpEnd_P;
  double b_gamma;
  double d;
  double delCLmax;
  double the_f;
  int HighAlpEnd_P;
  // ----------------- STALL ANGLES AND HIGH AOA ANGLES -----------------%
  //  AR_pre = [0.5,0.75,1,1.25,1.5,1.75,2, 3,4,6];
  //  alpStall_pre = [44,37,34,22,20,18,18, 18,17,17]*pi/180;
  //  alpLim_pre = [35,33,28,20,15,14,13, 13,13,13]*pi/180; % By Torres
  // alpLim_pre = [32,36,34,38,32,21,16,11,10,8]*pi/180; % Alp_Limit
  //  Alp_Limit
  HighAlpStart = coder::interp1(dv, dv1, 2.0920035851844676);
  // ----------------------- LIFT CURVE SLOPE ---------------------------%
  //  From McCormick
  //  -------- POTENTIAL LIFT AND VORTEX LIFT COEFFICIENTS --------------%
  // --------------------- FLAP DEFLECTION ------------------------------%
  d = cf / c;
  if (d == 1.0) {
    the_f = 0.0;
    LowAlpEnd_P = 0.27122385289355255;
    LowAlpEnd_N = -0.27122385289355255;
    LowAlpStart_P = 2.8703688006962405;
    LowAlpStart_N = -2.8703688006962405;
    b_HighAlpEnd_P = 3.1415926535897931 - HighAlpStart;
    HighAlpStart_N = -HighAlpStart;
    HighAlpEnd_N = -(3.1415926535897931 - HighAlpStart);
    alp += def;
    b_gamma = 0.0;
  } else {
    the_f = std::acos(2.0 * cf / c - 1.0);
    // eta_pre = [0.81*0.7 0.8*0.7 0.685*0.7 0.535 0.455 0.405 0.37 0.345];
    for (HighAlpEnd_P = 0; HighAlpEnd_P < 8; HighAlpEnd_P++) {
      b_dv[HighAlpEnd_P] =
          0.17453292519943295 * static_cast<double>(HighAlpEnd_P);
    }
    HighAlpEnd_N = std::abs(def);
    the_f = 2.59283849675502 *
            (1.0 - (the_f - std::sin(the_f)) / 3.1415926535897931) *
            coder::c_interp1(b_dv, dv2, HighAlpEnd_N) * def;
    delCLmax = coder::b_interp1(dv3, dv4, d) * the_f;
    //  For Linear CL vs Alp curve
    //  alp0_eff1 = alp0 - delCL/CLAlp;
    //  CLmaxP1 = CLAlp*(alpStallP - alp0) + delCLmax;
    //  CLmaxN1 = CLAlp*(alpStallN - alp0) + delCLmax;
    //  alpStallP_eff1 = alp0_eff1 + CLmaxP1/CLAlp;
    //  alpStallN_eff1 = alp0_eff1 + CLmaxN1/CLAlp;
    //  For Non-Linear CL vs Alp curve
    // fun = @(x) Kp*sin(0-x)*(cos(0-x))^2 + Kv*abs(sin(0-x))*sin(0-x)*cos(0-x)
    // - delCL;
    the_f = coder::fzero(the_f);
    // fun = @(x) Kp*sin(x-alp0_eff)*(cos(x-alp0_eff))^2 +
    // Kv*abs(sin(x-alp0_eff))*sin(x-alp0_eff)*cos(x-alp0_eff) - CLmaxP;
    LowAlpEnd_P = coder::fzero(the_f, delCLmax + 0.8620384028016963);
    // alpStallP_eff = fzero(@myfun_alpStall_eff,x0,[],Kp,Kv,alp0_eff,CLmaxP);
    // fun = @(x) Kp*sin(x-alp0_eff)*(cos(x-alp0_eff))^2 +
    // Kv*abs(sin(x-alp0_eff))*sin(x-alp0_eff)*cos(x-alp0_eff) - CLmaxN;
    LowAlpEnd_N = coder::b_fzero(the_f, delCLmax + -0.8620384028016963);
    // alpStallN_eff = fzero(@myfun_alpStall_eff,x0,[],Kp,Kv,alp0_eff,CLmaxN);
    LowAlpStart_P = LowAlpEnd_N + 3.1415926535897931;
    LowAlpStart_N = LowAlpEnd_P - 3.1415926535897931;
    //  High AoA parameters using equivalent flat plate method
    CT2 = c - cf;
    b_gamma = std::asin(std::sin(def) * cf /
                        std::sqrt((CT2 * CT2 + cf * cf) +
                                  2.0 * CT2 * cf * std::cos(HighAlpEnd_N)));
    b_HighAlpEnd_P = 3.1415926535897931 - HighAlpStart;
    HighAlpStart_N = -HighAlpStart;
    HighAlpEnd_N = -(3.1415926535897931 - HighAlpStart);
  }
  //     %% Calculate Aerodynamic coefficients based on regime
  if ((alp >= -3.1415926535897931) && (alp <= LowAlpStart_N)) {
    double CL_tmp;
    double a_eff;
    double b_CL_tmp;
    a_eff = (alp - the_f) + 3.1415926535897931;
    CL_tmp = std::sin(a_eff);
    b_CL_tmp = std::abs(CL_tmp);
    the_f = std::cos(a_eff);
    *CL = 2.59283849675502 * CL_tmp * (the_f * the_f) +
          3.1415926535897931 * b_CL_tmp * CL_tmp * the_f;
    delCLmax = std::sin(a_eff);
    *CD = std::abs(2.59283849675502 * b_CL_tmp * CL_tmp * the_f +
                   3.1415926535897931 * rt_powd_snf(delCLmax, 3.0)) +
          0.02;
    *CN = -(2.59283849675502 * delCLmax * std::abs(std::cos(a_eff)) +
            3.1415926535897931 * std::abs(delCLmax) * delCLmax);
    // CM = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
    *CM = -(0.25 - 0.175 * (1.0 - std::abs(a_eff - 3.1415926535897931) /
                                      1.5707963267948966)) *
          *CN;
    //  ------------------------------------------------- %
  } else if ((alp > LowAlpStart_N) && (alp < HighAlpEnd_N)) {
    double CL_tmp;
    double a;
    double a_eff;
    double b_CL_tmp;
    double y;
    //  Linear End
    a_eff = (LowAlpStart_N - the_f) + 3.1415926535897931;
    a = std::cos(a_eff);
    b_HighAlpEnd_P = std::sin(a_eff);
    LowAlpEnd_N = std::cos(a_eff);
    the_f = std::sin(a_eff);
    delCLmax = std::abs(the_f);
    HighAlpStart_N = -(2.59283849675502 * the_f * LowAlpEnd_N +
                       3.1415926535897931 * delCLmax * the_f);
    // CM1 = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
    // CD1 = Cd0 + abs(CL1)*abs(tan(a_eff));
    // CN1 = -sqrt(CL1^2 + CD1^2);
    y = std::abs(a_eff - 3.1415926535897931);
    //  NonLinear End
    a_eff = HighAlpEnd_N + b_gamma;
    // - alp0_eff;
    //  Determining effective Cd90
    if ((d == 1.0) || (d == 0.0)) {
      b_gamma = Cd90;
    } else {
      if (HighAlpEnd_N < 0.0) {
        HighAlpEnd_P = -1;
      } else {
        HighAlpEnd_P = (HighAlpEnd_N > 0.0);
      }
      if (HighAlpEnd_P == 1) {
        b_CL_tmp = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (b_CL_tmp * b_CL_tmp) +
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      } else {
        b_CL_tmp = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (b_CL_tmp * b_CL_tmp) -
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      }
    }
    //  Hoerner Model with Lindenberg Correction
    b_CL_tmp = std::sin(std::abs(a_eff));
    b_gamma = -b_gamma *
              (1.0 / (0.44 * b_CL_tmp + 0.56) - 0.40987876493629172) * b_CL_tmp;
    b_CL_tmp = std::cos(a_eff);
    CT2 = 0.01 * b_CL_tmp;
    CL_tmp = std::sin(a_eff);
    //  Interpolation
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    CN1[0] = HighAlpStart_N;
    CN1[1] = b_gamma;
    *CN = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    CN1[0] = 2.59283849675502 * b_HighAlpEnd_P * (a * a) +
             3.1415926535897931 * std::abs(b_HighAlpEnd_P) * b_HighAlpEnd_P *
                 LowAlpEnd_N;
    CN1[1] = b_gamma * b_CL_tmp - CT2 * CL_tmp;
    *CL = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    CN1[0] = std::abs(2.59283849675502 * delCLmax * the_f * LowAlpEnd_N +
                      3.1415926535897931 * rt_powd_snf(the_f, 3.0)) +
             0.02;
    CN1[1] = b_gamma * CL_tmp + CT2 * b_CL_tmp;
    *CD = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    CN1[0] = -(0.25 - 0.175 * (1.0 - y / 1.5707963267948966)) * HighAlpStart_N;
    CN1[1] = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) *
             b_gamma;
    *CM = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    //  ------------------------------------------------- %
  } else if ((alp >= HighAlpEnd_N) && (alp <= HighAlpStart_N)) {
    double CL_tmp;
    double a_eff;
    //  Equivalent flat plate
    a_eff = alp + b_gamma;
    // - alp0_eff;
    //  Determining effective Cd90
    if ((d == 1.0) || (d == 0.0)) {
      b_gamma = Cd90;
    } else {
      if (alp < 0.0) {
        HighAlpEnd_P = -1;
      } else {
        HighAlpEnd_P = (alp > 0.0);
      }
      if (HighAlpEnd_P == 1) {
        double a;
        a = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (a * a) +
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      } else {
        CT2 = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (CT2 * CT2) - 0.00361111 * CT2) + Cd90;
      }
    }
    //  Hoerner Model
    the_f = std::sin(std::abs(a_eff));
    *CN =
        -b_gamma * (1.0 / (0.44 * the_f + 0.56) - 0.40987876493629172) * the_f;
    the_f = std::cos(a_eff);
    delCLmax = 0.01 * the_f;
    CL_tmp = std::sin(a_eff);
    *CL = *CN * the_f - delCLmax * CL_tmp;
    *CD = *CN * CL_tmp + delCLmax * the_f;
    *CM = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) * *CN;
    //  ------------------------------------------------- %
  } else if ((alp > HighAlpStart_N) && (alp < LowAlpEnd_N)) {
    double CL_tmp;
    double CN1_tmp;
    double a;
    double a_eff;
    double b_CL_tmp;
    double y;
    //  Linear End
    a_eff = LowAlpEnd_N - the_f;
    LowAlpStart_N = std::sin(a_eff);
    the_f = std::abs(LowAlpStart_N);
    delCLmax = std::cos(a_eff);
    a = std::cos(a_eff);
    HighAlpEnd_N = rt_powd_snf(std::sin(a_eff), 3.0);
    CL_tmp = std::sin(a_eff);
    y = std::abs(std::sin(a_eff));
    b_HighAlpEnd_P = std::sin(a_eff);
    //  NonLinear End
    a_eff = HighAlpStart_N + b_gamma;
    // - alp0_eff;
    if ((d == 1.0) || (d == 0.0)) {
      b_gamma = Cd90;
    } else {
      if (HighAlpStart_N < 0.0) {
        HighAlpEnd_P = -1;
      } else {
        HighAlpEnd_P = (HighAlpStart_N > 0.0);
      }
      if (HighAlpEnd_P == 1) {
        b_CL_tmp = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (b_CL_tmp * b_CL_tmp) +
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      } else {
        b_CL_tmp = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (b_CL_tmp * b_CL_tmp) -
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      }
    }
    b_CL_tmp = std::sin(std::abs(a_eff));
    b_gamma = -b_gamma *
              (1.0 / (0.44 * b_CL_tmp + 0.56) - 0.40987876493629172) * b_CL_tmp;
    b_CL_tmp = std::cos(a_eff);
    CT2 = 0.01 * b_CL_tmp;
    //  Interpolation
    b_LowAlpStart_P[0] = LowAlpEnd_N;
    b_LowAlpStart_P[1] = HighAlpStart_N;
    CN1[0] = 2.59283849675502 * CL_tmp * delCLmax +
             3.1415926535897931 * y * b_HighAlpEnd_P;
    CN1[1] = b_gamma;
    *CN = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_N;
    b_LowAlpStart_P[1] = HighAlpStart_N;
    CN1[0] = 2.59283849675502 * LowAlpStart_N * (a * a) +
             3.1415926535897931 * the_f * LowAlpStart_N * delCLmax;
    CN1_tmp = std::sin(a_eff);
    CN1[1] = b_gamma * b_CL_tmp - CT2 * CN1_tmp;
    *CL = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_N;
    b_LowAlpStart_P[1] = HighAlpStart_N;
    CN1[0] = std::abs(2.59283849675502 * the_f * LowAlpStart_N * delCLmax +
                      3.1415926535897931 * HighAlpEnd_N) +
             0.02;
    CN1[1] = b_gamma * CN1_tmp + CT2 * std::cos(a_eff);
    *CD = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_N;
    b_LowAlpStart_P[1] = HighAlpStart_N;
    CN1[0] = -0.53407075111026481 * std::abs(b_HighAlpEnd_P) * b_HighAlpEnd_P;
    CN1[1] = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) *
             b_gamma;
    *CM = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    //  ------------------------------------------------- %
  } else if ((alp >= LowAlpEnd_N) && (alp <= LowAlpEnd_P)) {
    double a_eff;
    a_eff = alp - the_f;
    CT2 = std::cos(a_eff);
    *CL =
        2.59283849675502 * std::sin(a_eff) * (CT2 * CT2) +
        3.1415926535897931 * std::abs(std::sin(a_eff)) * std::sin(a_eff) * CT2;
    delCLmax = std::sin(a_eff);
    the_f = std::abs(delCLmax);
    *CD = (std::abs(2.59283849675502 * the_f * delCLmax * std::cos(a_eff)) +
           0.02) +
          3.1415926535897931 * std::abs(rt_powd_snf(delCLmax, 3.0));
    *CN = 2.59283849675502 * delCLmax * std::cos(a_eff) +
          3.1415926535897931 * the_f * delCLmax;
    *CM = -0.53407075111026481 * std::abs(std::sin(a_eff)) * std::sin(a_eff);
    //  ------------------------------------------------- %
  } else if ((alp > LowAlpEnd_P) && (alp < HighAlpStart)) {
    double CL_tmp;
    double CN1_tmp;
    double a;
    double a_eff;
    double b_CL_tmp;
    double y;
    //  Linear End
    a_eff = LowAlpEnd_P - the_f;
    a = std::cos(a_eff);
    CL_tmp = std::sin(a_eff);
    y = std::abs(std::sin(a_eff));
    b_HighAlpEnd_P = std::sin(a_eff);
    LowAlpEnd_N = std::cos(a_eff);
    delCLmax = std::sin(a_eff);
    the_f = std::cos(a_eff);
    //  NonLinear End
    a_eff = HighAlpStart + b_gamma;
    // - alp0_eff;
    if ((d == 1.0) || (d == 0.0)) {
      b_gamma = Cd90;
    } else {
      if (HighAlpStart < 0.0) {
        HighAlpEnd_P = -1;
      } else {
        HighAlpEnd_P = (HighAlpStart > 0.0);
      }
      if (HighAlpEnd_P == 1) {
        b_CL_tmp = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (b_CL_tmp * b_CL_tmp) +
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      } else {
        b_CL_tmp = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (b_CL_tmp * b_CL_tmp) -
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      }
    }
    b_CL_tmp = std::sin(std::abs(a_eff));
    b_gamma = b_gamma * (1.0 / (0.44 * b_CL_tmp + 0.56) - 0.40987876493629172) *
              b_CL_tmp;
    b_CL_tmp = std::cos(a_eff);
    CT2 = 0.01 * b_CL_tmp;
    //  Interpolation
    b_LowAlpStart_P[0] = LowAlpEnd_P;
    b_LowAlpStart_P[1] = HighAlpStart;
    CN1_tmp = std::abs(delCLmax);
    CN1[0] = 2.59283849675502 * delCLmax * the_f +
             3.1415926535897931 * CN1_tmp * delCLmax;
    CN1[1] = b_gamma;
    *CN = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_P;
    b_LowAlpStart_P[1] = HighAlpStart;
    CN1[0] = 2.59283849675502 * CL_tmp * (a * a) +
             3.1415926535897931 * y * b_HighAlpEnd_P * LowAlpEnd_N;
    CN1[1] = b_gamma * b_CL_tmp - CT2 * std::sin(a_eff);
    *CL = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_P;
    b_LowAlpStart_P[1] = HighAlpStart;
    CN1[0] = std::abs(2.59283849675502 * std::abs(b_HighAlpEnd_P) *
                          b_HighAlpEnd_P * LowAlpEnd_N +
                      3.1415926535897931 * rt_powd_snf(b_HighAlpEnd_P, 3.0)) +
             0.02;
    CN1[1] = b_gamma * std::sin(a_eff) + CT2 * b_CL_tmp;
    *CD = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpEnd_P;
    b_LowAlpStart_P[1] = HighAlpStart;
    CN1[0] = -0.53407075111026481 * CN1_tmp * delCLmax;
    CN1[1] = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) *
             b_gamma;
    *CM = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    //  ------------------------------------------------- %
  } else if ((alp >= HighAlpStart) && (alp <= b_HighAlpEnd_P)) {
    double CL_tmp;
    double a_eff;
    //  Equivalent flat plate
    a_eff = alp + b_gamma;
    // - alp0_eff;
    //  Determining effective Cd90
    if ((d == 1.0) || (d == 0.0)) {
      b_gamma = Cd90;
    } else {
      if (alp < 0.0) {
        HighAlpEnd_P = -1;
      } else {
        HighAlpEnd_P = (alp > 0.0);
      }
      if (HighAlpEnd_P == 1) {
        double a;
        a = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (a * a) +
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      } else {
        double a;
        a = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (a * a) -
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      }
    }
    //  Hoerner Model
    the_f = std::sin(std::abs(a_eff));
    *CN = b_gamma * (1.0 / (0.44 * the_f + 0.56) - 0.40987876493629172) * the_f;
    the_f = std::cos(a_eff);
    delCLmax = 0.01 * the_f;
    CL_tmp = std::sin(a_eff);
    *CL = *CN * the_f - delCLmax * CL_tmp;
    *CD = *CN * CL_tmp + delCLmax * the_f;
    *CM = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) * *CN;
    //  ------------------------------------------------- %
  } else if ((alp > b_HighAlpEnd_P) && (alp < LowAlpStart_P)) {
    double CL_tmp;
    double CN1_tmp;
    double a_eff;
    double b_CL_tmp;
    double y;
    //  Linear End
    a_eff = (LowAlpStart_P - the_f) - 3.1415926535897931;
    LowAlpStart_N = std::sin(a_eff);
    CT2 = std::cos(a_eff);
    y = std::abs(std::sin(a_eff));
    b_CL_tmp = std::sin(a_eff);
    LowAlpEnd_N = std::cos(a_eff);
    CN1_tmp = std::abs(b_CL_tmp);
    HighAlpStart_N = -(2.59283849675502 * b_CL_tmp * LowAlpEnd_N +
                       3.1415926535897931 * CN1_tmp * b_CL_tmp);
    HighAlpEnd_N = std::abs(a_eff - 3.1415926535897931);
    // CM1 = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
    //  NonLinear End
    a_eff = b_HighAlpEnd_P + b_gamma;
    // - alp0_eff;
    //  Determining effective Cd90
    if ((d == 1.0) || (d == 0.0)) {
      b_gamma = Cd90;
    } else {
      if (b_HighAlpEnd_P < 0.0) {
        HighAlpEnd_P = -1;
      } else {
        HighAlpEnd_P = (b_HighAlpEnd_P > 0.0);
      }
      if (HighAlpEnd_P == 1) {
        double a;
        a = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (a * a) +
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      } else {
        double a;
        a = def * 180.0 / 3.1415926535897931;
        b_gamma = (-1.2963E-5 * (a * a) -
                   0.00361111 * (def * 180.0 / 3.1415926535897931)) +
                  Cd90;
      }
    }
    //  Hoerner Model
    b_gamma = b_gamma *
              (1.0 / (0.44 * std::sin(std::abs(a_eff)) + 0.56) -
               0.40987876493629172) *
              std::sin(std::abs(a_eff));
    delCLmax = 0.01 * std::cos(a_eff);
    CL_tmp = std::sin(a_eff);
    the_f = std::cos(a_eff);
    //  Interpolation
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = b_HighAlpEnd_P;
    CN1[0] = HighAlpStart_N;
    CN1[1] = b_gamma;
    *CN = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = b_HighAlpEnd_P;
    CN1[0] = 2.59283849675502 * LowAlpStart_N * (CT2 * CT2) +
             3.1415926535897931 * y * LowAlpStart_N * CT2;
    CN1[1] = b_gamma * the_f - delCLmax * CL_tmp;
    *CL = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = b_HighAlpEnd_P;
    CN1[0] = std::abs(2.59283849675502 * CN1_tmp * b_CL_tmp * LowAlpEnd_N +
                      3.1415926535897931 * rt_powd_snf(b_CL_tmp, 3.0)) +
             0.02;
    CN1[1] = b_gamma * CL_tmp + delCLmax * the_f;
    *CD = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = b_HighAlpEnd_P;
    CN1[0] = -(0.25 - 0.175 * (1.0 - HighAlpEnd_N / 1.5707963267948966)) *
             HighAlpStart_N;
    CN1[1] = -(0.25 - 0.175 * (1.0 - std::abs(a_eff) / 1.5707963267948966)) *
             b_gamma;
    *CM = coder::d_interp1(b_LowAlpStart_P, CN1, alp);
    //  ------------------------------------------------- %
  } else {
    double CL_tmp;
    double a;
    double a_eff;
    double b_CL_tmp;
    a_eff = (alp - the_f) - 3.1415926535897931;
    a = std::cos(a_eff);
    CL_tmp = std::sin(a_eff);
    b_CL_tmp = std::cos(a_eff);
    the_f = std::abs(CL_tmp);
    delCLmax = 2.59283849675502 * CL_tmp;
    HighAlpEnd_N = 3.1415926535897931 * the_f * CL_tmp;
    *CL = delCLmax * (a * a) + HighAlpEnd_N * b_CL_tmp;
    *CD = std::abs(2.59283849675502 * the_f * CL_tmp * b_CL_tmp +
                   3.1415926535897931 * rt_powd_snf(CL_tmp, 3.0)) +
          0.02;
    *CN = -(delCLmax * b_CL_tmp + HighAlpEnd_N);
    *CM = -(0.25 - 0.175 * (1.0 - std::abs(a_eff - 3.1415926535897931) /
                                      1.5707963267948966)) *
          *CN;
    // CM = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
  }
}

//
// File trailer for flappedAirfoil_LAR_v1_FlatPlate.cpp
//
// [EOF]
//


//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: fzero.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 05-Aug-2021 13:41:18
//

// Include Files
#include "fzero.h"
#include "rt_nonfinite.h"
#include <cmath>

// Function Definitions
//
// Arguments    : double varargin_3
// Return Type  : double
//
namespace coder {
double b_fzero(double varargin_3)
{
  double b;
  if (0.0 - varargin_3 == 0.0) {
    b = 0.0;
  } else {
    double a;
    double dx;
    double fa;
    double fb;
    double q;
    double r;
    int exitg2;
    dx = 0.02;
    a = 0.0;
    fa = 0.0 - varargin_3;
    b = 0.0;
    fb = 0.0 - varargin_3;
    do {
      exitg2 = 0;
      if ((fa > 0.0) == (fb > 0.0)) {
        dx *= 1.4142135623730951;
        a = 0.0 - dx;
        r = std::sin(0.0 - (0.0 - dx));
        q = std::cos(0.0 - (0.0 - dx));
        fa = (1.7436202582132567 * r * (q * q) +
              3.1415926535897931 * std::abs(r) * r * q) -
             varargin_3;
        if (std::isinf(fa) || std::isnan(fa)) {
          b = rtNaN;
          exitg2 = 1;
        } else if (std::isinf(0.0 - dx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if ((fa > 0.0) != (fb > 0.0)) {
          exitg2 = 2;
        } else {
          b = dx;
          r = std::sin(0.0 - dx);
          q = std::cos(0.0 - dx);
          fb = (1.7436202582132567 * r * (q * q) +
                3.1415926535897931 * std::abs(r) * r * q) -
               varargin_3;
          if (std::isinf(fb) || std::isnan(fb)) {
            b = rtNaN;
            exitg2 = 1;
          } else if (std::isinf(dx)) {
            b = rtNaN;
            exitg2 = 1;
          }
        }
      } else {
        exitg2 = 2;
      }
    } while (exitg2 == 0);
    if (exitg2 != 1) {
      double c;
      double d;
      double e;
      double fc;
      boolean_T exitg1;
      fc = fb;
      c = b;
      e = 0.0;
      d = 0.0;
      exitg1 = false;
      while ((!exitg1) && ((fb != 0.0) && (a != b))) {
        double m;
        double toler;
        if ((fb > 0.0) == (fc > 0.0)) {
          c = a;
          fc = fa;
          d = b - a;
          e = d;
        }
        if (std::abs(fc) < std::abs(fb)) {
          a = b;
          b = c;
          c = a;
          fa = fb;
          fb = fc;
          fc = fa;
        }
        m = 0.5 * (c - b);
        toler = 4.4408920985006262E-16 * std::fmax(std::abs(b), 1.0);
        if ((std::abs(m) <= toler) || (fb == 0.0)) {
          exitg1 = true;
        } else {
          if ((std::abs(e) < toler) || (std::abs(fa) <= std::abs(fb))) {
            d = m;
            e = m;
          } else {
            double s;
            s = fb / fa;
            if (a == c) {
              dx = 2.0 * m * s;
              q = 1.0 - s;
            } else {
              q = fa / fc;
              r = fb / fc;
              dx = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
              q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (dx > 0.0) {
              q = -q;
            } else {
              dx = -dx;
            }
            if ((2.0 * dx < 3.0 * m * q - std::abs(toler * q)) &&
                (dx < std::abs(0.5 * e * q))) {
              e = d;
              d = dx / q;
            } else {
              d = m;
              e = m;
            }
          }
          a = b;
          fa = fb;
          if (std::abs(d) > toler) {
            b += d;
          } else if (b > c) {
            b -= toler;
          } else {
            b += toler;
          }
          r = std::sin(0.0 - b);
          dx = std::cos(0.0 - b);
          fb = (1.7436202582132567 * r * (dx * dx) +
                3.1415926535897931 * std::abs(r) * r * std::cos(0.0 - b)) -
               varargin_3;
        }
      }
    }
  }
  return b;
}

//
// Arguments    : double varargin_3
//                double varargin_4
// Return Type  : double
//
double b_fzero(double varargin_3, double varargin_4)
{
  double a_eff;
  double b;
  double fx;
  double r;
  r = std::sin(-0.27122385289355255 - varargin_3);
  a_eff = std::cos(-0.27122385289355255 - varargin_3);
  fx = (2.59283849675502 * r * (a_eff * a_eff) +
        3.1415926535897931 * std::abs(r) * r * a_eff) -
       varargin_4;
  if (fx == 0.0) {
    b = -0.27122385289355255;
  } else {
    double a;
    double dx;
    double fb;
    int exitg2;
    dx = -0.0054244770578710513;
    a = -0.27122385289355255;
    b = -0.27122385289355255;
    fb = fx;
    do {
      exitg2 = 0;
      if ((fx > 0.0) == (fb > 0.0)) {
        dx *= 1.4142135623730951;
        a = -0.27122385289355255 - dx;
        a_eff = (-0.27122385289355255 - dx) - varargin_3;
        r = std::sin(a_eff);
        a_eff = std::cos(a_eff);
        fx = (2.59283849675502 * r * (a_eff * a_eff) +
              3.1415926535897931 * std::abs(r) * r * a_eff) -
             varargin_4;
        if (std::isinf(fx) || std::isnan(fx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if (std::isinf(-0.27122385289355255 - dx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if ((fx > 0.0) != (fb > 0.0)) {
          exitg2 = 2;
        } else {
          b = dx + -0.27122385289355255;
          a_eff = (dx + -0.27122385289355255) - varargin_3;
          r = std::sin(a_eff);
          a_eff = std::cos(a_eff);
          fb = (2.59283849675502 * r * (a_eff * a_eff) +
                3.1415926535897931 * std::abs(r) * r * a_eff) -
               varargin_4;
          if (std::isinf(fb) || std::isnan(fb)) {
            b = rtNaN;
            exitg2 = 1;
          } else if (std::isinf(dx + -0.27122385289355255)) {
            b = rtNaN;
            exitg2 = 1;
          }
        }
      } else {
        exitg2 = 2;
      }
    } while (exitg2 == 0);
    if (exitg2 != 1) {
      double c;
      double d;
      double e;
      double fc;
      boolean_T exitg1;
      fc = fb;
      c = b;
      e = 0.0;
      d = 0.0;
      exitg1 = false;
      while ((!exitg1) && ((fb != 0.0) && (a != b))) {
        double m;
        double toler;
        if ((fb > 0.0) == (fc > 0.0)) {
          c = a;
          fc = fx;
          d = b - a;
          e = d;
        }
        if (std::abs(fc) < std::abs(fb)) {
          a = b;
          b = c;
          c = a;
          fx = fb;
          fb = fc;
          fc = fx;
        }
        m = 0.5 * (c - b);
        toler = 4.4408920985006262E-16 * std::fmax(std::abs(b), 1.0);
        if ((std::abs(m) <= toler) || (fb == 0.0)) {
          exitg1 = true;
        } else {
          if ((std::abs(e) < toler) || (std::abs(fx) <= std::abs(fb))) {
            d = m;
            e = m;
          } else {
            double s;
            s = fb / fx;
            if (a == c) {
              dx = 2.0 * m * s;
              a_eff = 1.0 - s;
            } else {
              a_eff = fx / fc;
              r = fb / fc;
              dx = s * (2.0 * m * a_eff * (a_eff - r) - (b - a) * (r - 1.0));
              a_eff = (a_eff - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (dx > 0.0) {
              a_eff = -a_eff;
            } else {
              dx = -dx;
            }
            if ((2.0 * dx < 3.0 * m * a_eff - std::abs(toler * a_eff)) &&
                (dx < std::abs(0.5 * e * a_eff))) {
              e = d;
              d = dx / a_eff;
            } else {
              d = m;
              e = m;
            }
          }
          a = b;
          fx = fb;
          if (std::abs(d) > toler) {
            b += d;
          } else if (b > c) {
            b -= toler;
          } else {
            b += toler;
          }
          a_eff = b - varargin_3;
          r = std::sin(a_eff);
          a_eff = std::cos(a_eff);
          fb = (2.59283849675502 * r * (a_eff * a_eff) +
                3.1415926535897931 * std::abs(r) * r * a_eff) -
               varargin_4;
        }
      }
    }
  }
  return b;
}

//
// Arguments    : double varargin_3
//                double varargin_4
// Return Type  : double
//
double c_fzero(double varargin_3, double varargin_4)
{
  double a_eff;
  double b;
  double fx;
  double r;
  r = std::sin(0.37260508566838529 - varargin_3);
  a_eff = std::cos(0.37260508566838529 - varargin_3);
  fx = (1.7436202582132567 * r * (a_eff * a_eff) +
        3.1415926535897931 * std::abs(r) * r * a_eff) -
       varargin_4;
  if (fx == 0.0) {
    b = 0.37260508566838529;
  } else {
    double a;
    double dx;
    double fb;
    int exitg2;
    dx = 0.0074521017133677061;
    a = 0.37260508566838529;
    b = 0.37260508566838529;
    fb = fx;
    do {
      exitg2 = 0;
      if ((fx > 0.0) == (fb > 0.0)) {
        dx *= 1.4142135623730951;
        a = 0.37260508566838529 - dx;
        a_eff = (0.37260508566838529 - dx) - varargin_3;
        r = std::sin(a_eff);
        a_eff = std::cos(a_eff);
        fx = (1.7436202582132567 * r * (a_eff * a_eff) +
              3.1415926535897931 * std::abs(r) * r * a_eff) -
             varargin_4;
        if (std::isinf(fx) || std::isnan(fx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if (std::isinf(0.37260508566838529 - dx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if ((fx > 0.0) != (fb > 0.0)) {
          exitg2 = 2;
        } else {
          b = dx + 0.37260508566838529;
          a_eff = (dx + 0.37260508566838529) - varargin_3;
          r = std::sin(a_eff);
          a_eff = std::cos(a_eff);
          fb = (1.7436202582132567 * r * (a_eff * a_eff) +
                3.1415926535897931 * std::abs(r) * r * a_eff) -
               varargin_4;
          if (std::isinf(fb) || std::isnan(fb)) {
            b = rtNaN;
            exitg2 = 1;
          } else if (std::isinf(dx + 0.37260508566838529)) {
            b = rtNaN;
            exitg2 = 1;
          }
        }
      } else {
        exitg2 = 2;
      }
    } while (exitg2 == 0);
    if (exitg2 != 1) {
      double c;
      double d;
      double e;
      double fc;
      boolean_T exitg1;
      fc = fb;
      c = b;
      e = 0.0;
      d = 0.0;
      exitg1 = false;
      while ((!exitg1) && ((fb != 0.0) && (a != b))) {
        double m;
        double toler;
        if ((fb > 0.0) == (fc > 0.0)) {
          c = a;
          fc = fx;
          d = b - a;
          e = d;
        }
        if (std::abs(fc) < std::abs(fb)) {
          a = b;
          b = c;
          c = a;
          fx = fb;
          fb = fc;
          fc = fx;
        }
        m = 0.5 * (c - b);
        toler = 4.4408920985006262E-16 * std::fmax(std::abs(b), 1.0);
        if ((std::abs(m) <= toler) || (fb == 0.0)) {
          exitg1 = true;
        } else {
          if ((std::abs(e) < toler) || (std::abs(fx) <= std::abs(fb))) {
            d = m;
            e = m;
          } else {
            double s;
            s = fb / fx;
            if (a == c) {
              dx = 2.0 * m * s;
              a_eff = 1.0 - s;
            } else {
              a_eff = fx / fc;
              r = fb / fc;
              dx = s * (2.0 * m * a_eff * (a_eff - r) - (b - a) * (r - 1.0));
              a_eff = (a_eff - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (dx > 0.0) {
              a_eff = -a_eff;
            } else {
              dx = -dx;
            }
            if ((2.0 * dx < 3.0 * m * a_eff - std::abs(toler * a_eff)) &&
                (dx < std::abs(0.5 * e * a_eff))) {
              e = d;
              d = dx / a_eff;
            } else {
              d = m;
              e = m;
            }
          }
          a = b;
          fx = fb;
          if (std::abs(d) > toler) {
            b += d;
          } else if (b > c) {
            b -= toler;
          } else {
            b += toler;
          }
          a_eff = b - varargin_3;
          r = std::sin(a_eff);
          a_eff = std::cos(a_eff);
          fb = (1.7436202582132567 * r * (a_eff * a_eff) +
                3.1415926535897931 * std::abs(r) * r * a_eff) -
               varargin_4;
        }
      }
    }
  }
  return b;
}

//
// Arguments    : double varargin_3
// Return Type  : double
//
double c_fzero(double varargin_3)
{
  double b;
  if (0.0 - varargin_3 == 0.0) {
    b = 0.0;
  } else {
    double a;
    double dx;
    double fa;
    double fb;
    double q;
    double r;
    int exitg2;
    dx = 0.02;
    a = 0.0;
    fa = 0.0 - varargin_3;
    b = 0.0;
    fb = 0.0 - varargin_3;
    do {
      exitg2 = 0;
      if ((fa > 0.0) == (fb > 0.0)) {
        dx *= 1.4142135623730951;
        a = 0.0 - dx;
        r = std::sin(0.0 - (0.0 - dx));
        q = std::cos(0.0 - (0.0 - dx));
        fa = (1.8084579108343679 * r * (q * q) +
              3.1415926535897931 * std::abs(r) * r * q) -
             varargin_3;
        if (std::isinf(fa) || std::isnan(fa)) {
          b = rtNaN;
          exitg2 = 1;
        } else if (std::isinf(0.0 - dx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if ((fa > 0.0) != (fb > 0.0)) {
          exitg2 = 2;
        } else {
          b = dx;
          r = std::sin(0.0 - dx);
          q = std::cos(0.0 - dx);
          fb = (1.8084579108343679 * r * (q * q) +
                3.1415926535897931 * std::abs(r) * r * q) -
               varargin_3;
          if (std::isinf(fb) || std::isnan(fb)) {
            b = rtNaN;
            exitg2 = 1;
          } else if (std::isinf(dx)) {
            b = rtNaN;
            exitg2 = 1;
          }
        }
      } else {
        exitg2 = 2;
      }
    } while (exitg2 == 0);
    if (exitg2 != 1) {
      double c;
      double d;
      double e;
      double fc;
      boolean_T exitg1;
      fc = fb;
      c = b;
      e = 0.0;
      d = 0.0;
      exitg1 = false;
      while ((!exitg1) && ((fb != 0.0) && (a != b))) {
        double m;
        double toler;
        if ((fb > 0.0) == (fc > 0.0)) {
          c = a;
          fc = fa;
          d = b - a;
          e = d;
        }
        if (std::abs(fc) < std::abs(fb)) {
          a = b;
          b = c;
          c = a;
          fa = fb;
          fb = fc;
          fc = fa;
        }
        m = 0.5 * (c - b);
        toler = 4.4408920985006262E-16 * std::fmax(std::abs(b), 1.0);
        if ((std::abs(m) <= toler) || (fb == 0.0)) {
          exitg1 = true;
        } else {
          if ((std::abs(e) < toler) || (std::abs(fa) <= std::abs(fb))) {
            d = m;
            e = m;
          } else {
            double s;
            s = fb / fa;
            if (a == c) {
              dx = 2.0 * m * s;
              q = 1.0 - s;
            } else {
              q = fa / fc;
              r = fb / fc;
              dx = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
              q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (dx > 0.0) {
              q = -q;
            } else {
              dx = -dx;
            }
            if ((2.0 * dx < 3.0 * m * q - std::abs(toler * q)) &&
                (dx < std::abs(0.5 * e * q))) {
              e = d;
              d = dx / q;
            } else {
              d = m;
              e = m;
            }
          }
          a = b;
          fa = fb;
          if (std::abs(d) > toler) {
            b += d;
          } else if (b > c) {
            b -= toler;
          } else {
            b += toler;
          }
          r = std::sin(0.0 - b);
          dx = std::cos(0.0 - b);
          fb = (1.8084579108343679 * r * (dx * dx) +
                3.1415926535897931 * std::abs(r) * r * std::cos(0.0 - b)) -
               varargin_3;
        }
      }
    }
  }
  return b;
}

//
// Arguments    : double varargin_3
//                double varargin_4
// Return Type  : double
//
double d_fzero(double varargin_3, double varargin_4)
{
  double a_eff;
  double b;
  double fx;
  double r;
  r = std::sin(-0.37260508566838529 - varargin_3);
  a_eff = std::cos(-0.37260508566838529 - varargin_3);
  fx = (1.7436202582132567 * r * (a_eff * a_eff) +
        3.1415926535897931 * std::abs(r) * r * a_eff) -
       varargin_4;
  if (fx == 0.0) {
    b = -0.37260508566838529;
  } else {
    double a;
    double dx;
    double fb;
    int exitg2;
    dx = -0.0074521017133677061;
    a = -0.37260508566838529;
    b = -0.37260508566838529;
    fb = fx;
    do {
      exitg2 = 0;
      if ((fx > 0.0) == (fb > 0.0)) {
        dx *= 1.4142135623730951;
        a = -0.37260508566838529 - dx;
        a_eff = (-0.37260508566838529 - dx) - varargin_3;
        r = std::sin(a_eff);
        a_eff = std::cos(a_eff);
        fx = (1.7436202582132567 * r * (a_eff * a_eff) +
              3.1415926535897931 * std::abs(r) * r * a_eff) -
             varargin_4;
        if (std::isinf(fx) || std::isnan(fx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if (std::isinf(-0.37260508566838529 - dx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if ((fx > 0.0) != (fb > 0.0)) {
          exitg2 = 2;
        } else {
          b = dx + -0.37260508566838529;
          a_eff = (dx + -0.37260508566838529) - varargin_3;
          r = std::sin(a_eff);
          a_eff = std::cos(a_eff);
          fb = (1.7436202582132567 * r * (a_eff * a_eff) +
                3.1415926535897931 * std::abs(r) * r * a_eff) -
               varargin_4;
          if (std::isinf(fb) || std::isnan(fb)) {
            b = rtNaN;
            exitg2 = 1;
          } else if (std::isinf(dx + -0.37260508566838529)) {
            b = rtNaN;
            exitg2 = 1;
          }
        }
      } else {
        exitg2 = 2;
      }
    } while (exitg2 == 0);
    if (exitg2 != 1) {
      double c;
      double d;
      double e;
      double fc;
      boolean_T exitg1;
      fc = fb;
      c = b;
      e = 0.0;
      d = 0.0;
      exitg1 = false;
      while ((!exitg1) && ((fb != 0.0) && (a != b))) {
        double m;
        double toler;
        if ((fb > 0.0) == (fc > 0.0)) {
          c = a;
          fc = fx;
          d = b - a;
          e = d;
        }
        if (std::abs(fc) < std::abs(fb)) {
          a = b;
          b = c;
          c = a;
          fx = fb;
          fb = fc;
          fc = fx;
        }
        m = 0.5 * (c - b);
        toler = 4.4408920985006262E-16 * std::fmax(std::abs(b), 1.0);
        if ((std::abs(m) <= toler) || (fb == 0.0)) {
          exitg1 = true;
        } else {
          if ((std::abs(e) < toler) || (std::abs(fx) <= std::abs(fb))) {
            d = m;
            e = m;
          } else {
            double s;
            s = fb / fx;
            if (a == c) {
              dx = 2.0 * m * s;
              a_eff = 1.0 - s;
            } else {
              a_eff = fx / fc;
              r = fb / fc;
              dx = s * (2.0 * m * a_eff * (a_eff - r) - (b - a) * (r - 1.0));
              a_eff = (a_eff - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (dx > 0.0) {
              a_eff = -a_eff;
            } else {
              dx = -dx;
            }
            if ((2.0 * dx < 3.0 * m * a_eff - std::abs(toler * a_eff)) &&
                (dx < std::abs(0.5 * e * a_eff))) {
              e = d;
              d = dx / a_eff;
            } else {
              d = m;
              e = m;
            }
          }
          a = b;
          fx = fb;
          if (std::abs(d) > toler) {
            b += d;
          } else if (b > c) {
            b -= toler;
          } else {
            b += toler;
          }
          a_eff = b - varargin_3;
          r = std::sin(a_eff);
          a_eff = std::cos(a_eff);
          fb = (1.7436202582132567 * r * (a_eff * a_eff) +
                3.1415926535897931 * std::abs(r) * r * a_eff) -
               varargin_4;
        }
      }
    }
  }
  return b;
}

//
// Arguments    : double varargin_3
// Return Type  : double
//
double d_fzero(double varargin_3)
{
  double b;
  if (0.0 - varargin_3 == 0.0) {
    b = 0.0;
  } else {
    double a;
    double dx;
    double fa;
    double fb;
    double q;
    double r;
    int exitg2;
    dx = 0.02;
    a = 0.0;
    fa = 0.0 - varargin_3;
    b = 0.0;
    fb = 0.0 - varargin_3;
    do {
      exitg2 = 0;
      if ((fa > 0.0) == (fb > 0.0)) {
        dx *= 1.4142135623730951;
        a = 0.0 - dx;
        r = std::sin(0.0 - (0.0 - dx));
        q = std::cos(0.0 - (0.0 - dx));
        fa = (0.35596863975551224 * r * (q * q) +
              3.1415926535897931 * std::abs(r) * r * q) -
             varargin_3;
        if (std::isinf(fa) || std::isnan(fa)) {
          b = rtNaN;
          exitg2 = 1;
        } else if (std::isinf(0.0 - dx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if ((fa > 0.0) != (fb > 0.0)) {
          exitg2 = 2;
        } else {
          b = dx;
          r = std::sin(0.0 - dx);
          q = std::cos(0.0 - dx);
          fb = (0.35596863975551224 * r * (q * q) +
                3.1415926535897931 * std::abs(r) * r * q) -
               varargin_3;
          if (std::isinf(fb) || std::isnan(fb)) {
            b = rtNaN;
            exitg2 = 1;
          } else if (std::isinf(dx)) {
            b = rtNaN;
            exitg2 = 1;
          }
        }
      } else {
        exitg2 = 2;
      }
    } while (exitg2 == 0);
    if (exitg2 != 1) {
      double c;
      double d;
      double e;
      double fc;
      boolean_T exitg1;
      fc = fb;
      c = b;
      e = 0.0;
      d = 0.0;
      exitg1 = false;
      while ((!exitg1) && ((fb != 0.0) && (a != b))) {
        double m;
        double toler;
        if ((fb > 0.0) == (fc > 0.0)) {
          c = a;
          fc = fa;
          d = b - a;
          e = d;
        }
        if (std::abs(fc) < std::abs(fb)) {
          a = b;
          b = c;
          c = a;
          fa = fb;
          fb = fc;
          fc = fa;
        }
        m = 0.5 * (c - b);
        toler = 4.4408920985006262E-16 * std::fmax(std::abs(b), 1.0);
        if ((std::abs(m) <= toler) || (fb == 0.0)) {
          exitg1 = true;
        } else {
          if ((std::abs(e) < toler) || (std::abs(fa) <= std::abs(fb))) {
            d = m;
            e = m;
          } else {
            double s;
            s = fb / fa;
            if (a == c) {
              dx = 2.0 * m * s;
              q = 1.0 - s;
            } else {
              q = fa / fc;
              r = fb / fc;
              dx = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
              q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (dx > 0.0) {
              q = -q;
            } else {
              dx = -dx;
            }
            if ((2.0 * dx < 3.0 * m * q - std::abs(toler * q)) &&
                (dx < std::abs(0.5 * e * q))) {
              e = d;
              d = dx / q;
            } else {
              d = m;
              e = m;
            }
          }
          a = b;
          fa = fb;
          if (std::abs(d) > toler) {
            b += d;
          } else if (b > c) {
            b -= toler;
          } else {
            b += toler;
          }
          r = std::sin(0.0 - b);
          dx = std::cos(0.0 - b);
          fb = (0.35596863975551224 * r * (dx * dx) +
                3.1415926535897931 * std::abs(r) * r * std::cos(0.0 - b)) -
               varargin_3;
        }
      }
    }
  }
  return b;
}

//
// Arguments    : double varargin_3
//                double varargin_4
// Return Type  : double
//
double e_fzero(double varargin_3, double varargin_4)
{
  double a_eff;
  double b;
  double fx;
  double r;
  r = std::sin(0.36088611410052907 - varargin_3);
  a_eff = std::cos(0.36088611410052907 - varargin_3);
  fx = (1.8084579108343679 * r * (a_eff * a_eff) +
        3.1415926535897931 * std::abs(r) * r * a_eff) -
       varargin_4;
  if (fx == 0.0) {
    b = 0.36088611410052907;
  } else {
    double a;
    double dx;
    double fb;
    int exitg2;
    dx = 0.0072177222820105809;
    a = 0.36088611410052907;
    b = 0.36088611410052907;
    fb = fx;
    do {
      exitg2 = 0;
      if ((fx > 0.0) == (fb > 0.0)) {
        dx *= 1.4142135623730951;
        a = 0.36088611410052907 - dx;
        a_eff = (0.36088611410052907 - dx) - varargin_3;
        r = std::sin(a_eff);
        a_eff = std::cos(a_eff);
        fx = (1.8084579108343679 * r * (a_eff * a_eff) +
              3.1415926535897931 * std::abs(r) * r * a_eff) -
             varargin_4;
        if (std::isinf(fx) || std::isnan(fx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if (std::isinf(0.36088611410052907 - dx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if ((fx > 0.0) != (fb > 0.0)) {
          exitg2 = 2;
        } else {
          b = dx + 0.36088611410052907;
          a_eff = (dx + 0.36088611410052907) - varargin_3;
          r = std::sin(a_eff);
          a_eff = std::cos(a_eff);
          fb = (1.8084579108343679 * r * (a_eff * a_eff) +
                3.1415926535897931 * std::abs(r) * r * a_eff) -
               varargin_4;
          if (std::isinf(fb) || std::isnan(fb)) {
            b = rtNaN;
            exitg2 = 1;
          } else if (std::isinf(dx + 0.36088611410052907)) {
            b = rtNaN;
            exitg2 = 1;
          }
        }
      } else {
        exitg2 = 2;
      }
    } while (exitg2 == 0);
    if (exitg2 != 1) {
      double c;
      double d;
      double e;
      double fc;
      boolean_T exitg1;
      fc = fb;
      c = b;
      e = 0.0;
      d = 0.0;
      exitg1 = false;
      while ((!exitg1) && ((fb != 0.0) && (a != b))) {
        double m;
        double toler;
        if ((fb > 0.0) == (fc > 0.0)) {
          c = a;
          fc = fx;
          d = b - a;
          e = d;
        }
        if (std::abs(fc) < std::abs(fb)) {
          a = b;
          b = c;
          c = a;
          fx = fb;
          fb = fc;
          fc = fx;
        }
        m = 0.5 * (c - b);
        toler = 4.4408920985006262E-16 * std::fmax(std::abs(b), 1.0);
        if ((std::abs(m) <= toler) || (fb == 0.0)) {
          exitg1 = true;
        } else {
          if ((std::abs(e) < toler) || (std::abs(fx) <= std::abs(fb))) {
            d = m;
            e = m;
          } else {
            double s;
            s = fb / fx;
            if (a == c) {
              dx = 2.0 * m * s;
              a_eff = 1.0 - s;
            } else {
              a_eff = fx / fc;
              r = fb / fc;
              dx = s * (2.0 * m * a_eff * (a_eff - r) - (b - a) * (r - 1.0));
              a_eff = (a_eff - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (dx > 0.0) {
              a_eff = -a_eff;
            } else {
              dx = -dx;
            }
            if ((2.0 * dx < 3.0 * m * a_eff - std::abs(toler * a_eff)) &&
                (dx < std::abs(0.5 * e * a_eff))) {
              e = d;
              d = dx / a_eff;
            } else {
              d = m;
              e = m;
            }
          }
          a = b;
          fx = fb;
          if (std::abs(d) > toler) {
            b += d;
          } else if (b > c) {
            b -= toler;
          } else {
            b += toler;
          }
          a_eff = b - varargin_3;
          r = std::sin(a_eff);
          a_eff = std::cos(a_eff);
          fb = (1.8084579108343679 * r * (a_eff * a_eff) +
                3.1415926535897931 * std::abs(r) * r * a_eff) -
               varargin_4;
        }
      }
    }
  }
  return b;
}

//
// Arguments    : double varargin_3
//                double varargin_4
// Return Type  : double
//
double f_fzero(double varargin_3, double varargin_4)
{
  double a_eff;
  double b;
  double fx;
  double r;
  r = std::sin(-0.36088611410052907 - varargin_3);
  a_eff = std::cos(-0.36088611410052907 - varargin_3);
  fx = (1.8084579108343679 * r * (a_eff * a_eff) +
        3.1415926535897931 * std::abs(r) * r * a_eff) -
       varargin_4;
  if (fx == 0.0) {
    b = -0.36088611410052907;
  } else {
    double a;
    double dx;
    double fb;
    int exitg2;
    dx = -0.0072177222820105809;
    a = -0.36088611410052907;
    b = -0.36088611410052907;
    fb = fx;
    do {
      exitg2 = 0;
      if ((fx > 0.0) == (fb > 0.0)) {
        dx *= 1.4142135623730951;
        a = -0.36088611410052907 - dx;
        a_eff = (-0.36088611410052907 - dx) - varargin_3;
        r = std::sin(a_eff);
        a_eff = std::cos(a_eff);
        fx = (1.8084579108343679 * r * (a_eff * a_eff) +
              3.1415926535897931 * std::abs(r) * r * a_eff) -
             varargin_4;
        if (std::isinf(fx) || std::isnan(fx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if (std::isinf(-0.36088611410052907 - dx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if ((fx > 0.0) != (fb > 0.0)) {
          exitg2 = 2;
        } else {
          b = dx + -0.36088611410052907;
          a_eff = (dx + -0.36088611410052907) - varargin_3;
          r = std::sin(a_eff);
          a_eff = std::cos(a_eff);
          fb = (1.8084579108343679 * r * (a_eff * a_eff) +
                3.1415926535897931 * std::abs(r) * r * a_eff) -
               varargin_4;
          if (std::isinf(fb) || std::isnan(fb)) {
            b = rtNaN;
            exitg2 = 1;
          } else if (std::isinf(dx + -0.36088611410052907)) {
            b = rtNaN;
            exitg2 = 1;
          }
        }
      } else {
        exitg2 = 2;
      }
    } while (exitg2 == 0);
    if (exitg2 != 1) {
      double c;
      double d;
      double e;
      double fc;
      boolean_T exitg1;
      fc = fb;
      c = b;
      e = 0.0;
      d = 0.0;
      exitg1 = false;
      while ((!exitg1) && ((fb != 0.0) && (a != b))) {
        double m;
        double toler;
        if ((fb > 0.0) == (fc > 0.0)) {
          c = a;
          fc = fx;
          d = b - a;
          e = d;
        }
        if (std::abs(fc) < std::abs(fb)) {
          a = b;
          b = c;
          c = a;
          fx = fb;
          fb = fc;
          fc = fx;
        }
        m = 0.5 * (c - b);
        toler = 4.4408920985006262E-16 * std::fmax(std::abs(b), 1.0);
        if ((std::abs(m) <= toler) || (fb == 0.0)) {
          exitg1 = true;
        } else {
          if ((std::abs(e) < toler) || (std::abs(fx) <= std::abs(fb))) {
            d = m;
            e = m;
          } else {
            double s;
            s = fb / fx;
            if (a == c) {
              dx = 2.0 * m * s;
              a_eff = 1.0 - s;
            } else {
              a_eff = fx / fc;
              r = fb / fc;
              dx = s * (2.0 * m * a_eff * (a_eff - r) - (b - a) * (r - 1.0));
              a_eff = (a_eff - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (dx > 0.0) {
              a_eff = -a_eff;
            } else {
              dx = -dx;
            }
            if ((2.0 * dx < 3.0 * m * a_eff - std::abs(toler * a_eff)) &&
                (dx < std::abs(0.5 * e * a_eff))) {
              e = d;
              d = dx / a_eff;
            } else {
              d = m;
              e = m;
            }
          }
          a = b;
          fx = fb;
          if (std::abs(d) > toler) {
            b += d;
          } else if (b > c) {
            b -= toler;
          } else {
            b += toler;
          }
          a_eff = b - varargin_3;
          r = std::sin(a_eff);
          a_eff = std::cos(a_eff);
          fb = (1.8084579108343679 * r * (a_eff * a_eff) +
                3.1415926535897931 * std::abs(r) * r * a_eff) -
               varargin_4;
        }
      }
    }
  }
  return b;
}

//
// Arguments    : double varargin_3
// Return Type  : double
//
double fzero(double varargin_3)
{
  double b;
  if (0.0 - varargin_3 == 0.0) {
    b = 0.0;
  } else {
    double a;
    double dx;
    double fa;
    double fb;
    double q;
    double r;
    int exitg2;
    dx = 0.02;
    a = 0.0;
    fa = 0.0 - varargin_3;
    b = 0.0;
    fb = 0.0 - varargin_3;
    do {
      exitg2 = 0;
      if ((fa > 0.0) == (fb > 0.0)) {
        dx *= 1.4142135623730951;
        a = 0.0 - dx;
        r = std::sin(0.0 - (0.0 - dx));
        q = std::cos(0.0 - (0.0 - dx));
        fa = (2.59283849675502 * r * (q * q) +
              3.1415926535897931 * std::abs(r) * r * q) -
             varargin_3;
        if (std::isinf(fa) || std::isnan(fa)) {
          b = rtNaN;
          exitg2 = 1;
        } else if (std::isinf(0.0 - dx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if ((fa > 0.0) != (fb > 0.0)) {
          exitg2 = 2;
        } else {
          b = dx;
          r = std::sin(0.0 - dx);
          q = std::cos(0.0 - dx);
          fb = (2.59283849675502 * r * (q * q) +
                3.1415926535897931 * std::abs(r) * r * q) -
               varargin_3;
          if (std::isinf(fb) || std::isnan(fb)) {
            b = rtNaN;
            exitg2 = 1;
          } else if (std::isinf(dx)) {
            b = rtNaN;
            exitg2 = 1;
          }
        }
      } else {
        exitg2 = 2;
      }
    } while (exitg2 == 0);
    if (exitg2 != 1) {
      double c;
      double d;
      double e;
      double fc;
      boolean_T exitg1;
      fc = fb;
      c = b;
      e = 0.0;
      d = 0.0;
      exitg1 = false;
      while ((!exitg1) && ((fb != 0.0) && (a != b))) {
        double m;
        double toler;
        if ((fb > 0.0) == (fc > 0.0)) {
          c = a;
          fc = fa;
          d = b - a;
          e = d;
        }
        if (std::abs(fc) < std::abs(fb)) {
          a = b;
          b = c;
          c = a;
          fa = fb;
          fb = fc;
          fc = fa;
        }
        m = 0.5 * (c - b);
        toler = 4.4408920985006262E-16 * std::fmax(std::abs(b), 1.0);
        if ((std::abs(m) <= toler) || (fb == 0.0)) {
          exitg1 = true;
        } else {
          if ((std::abs(e) < toler) || (std::abs(fa) <= std::abs(fb))) {
            d = m;
            e = m;
          } else {
            double s;
            s = fb / fa;
            if (a == c) {
              dx = 2.0 * m * s;
              q = 1.0 - s;
            } else {
              q = fa / fc;
              r = fb / fc;
              dx = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
              q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (dx > 0.0) {
              q = -q;
            } else {
              dx = -dx;
            }
            if ((2.0 * dx < 3.0 * m * q - std::abs(toler * q)) &&
                (dx < std::abs(0.5 * e * q))) {
              e = d;
              d = dx / q;
            } else {
              d = m;
              e = m;
            }
          }
          a = b;
          fa = fb;
          if (std::abs(d) > toler) {
            b += d;
          } else if (b > c) {
            b -= toler;
          } else {
            b += toler;
          }
          r = std::sin(0.0 - b);
          dx = std::cos(0.0 - b);
          fb = (2.59283849675502 * r * (dx * dx) +
                3.1415926535897931 * std::abs(r) * r * std::cos(0.0 - b)) -
               varargin_3;
        }
      }
    }
  }
  return b;
}

//
// Arguments    : double varargin_3
//                double varargin_4
// Return Type  : double
//
double fzero(double varargin_3, double varargin_4)
{
  double a_eff;
  double b;
  double fx;
  double r;
  r = std::sin(0.27122385289355255 - varargin_3);
  a_eff = std::cos(0.27122385289355255 - varargin_3);
  fx = (2.59283849675502 * r * (a_eff * a_eff) +
        3.1415926535897931 * std::abs(r) * r * a_eff) -
       varargin_4;
  if (fx == 0.0) {
    b = 0.27122385289355255;
  } else {
    double a;
    double dx;
    double fb;
    int exitg2;
    dx = 0.0054244770578710513;
    a = 0.27122385289355255;
    b = 0.27122385289355255;
    fb = fx;
    do {
      exitg2 = 0;
      if ((fx > 0.0) == (fb > 0.0)) {
        dx *= 1.4142135623730951;
        a = 0.27122385289355255 - dx;
        a_eff = (0.27122385289355255 - dx) - varargin_3;
        r = std::sin(a_eff);
        a_eff = std::cos(a_eff);
        fx = (2.59283849675502 * r * (a_eff * a_eff) +
              3.1415926535897931 * std::abs(r) * r * a_eff) -
             varargin_4;
        if (std::isinf(fx) || std::isnan(fx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if (std::isinf(0.27122385289355255 - dx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if ((fx > 0.0) != (fb > 0.0)) {
          exitg2 = 2;
        } else {
          b = dx + 0.27122385289355255;
          a_eff = (dx + 0.27122385289355255) - varargin_3;
          r = std::sin(a_eff);
          a_eff = std::cos(a_eff);
          fb = (2.59283849675502 * r * (a_eff * a_eff) +
                3.1415926535897931 * std::abs(r) * r * a_eff) -
               varargin_4;
          if (std::isinf(fb) || std::isnan(fb)) {
            b = rtNaN;
            exitg2 = 1;
          } else if (std::isinf(dx + 0.27122385289355255)) {
            b = rtNaN;
            exitg2 = 1;
          }
        }
      } else {
        exitg2 = 2;
      }
    } while (exitg2 == 0);
    if (exitg2 != 1) {
      double c;
      double d;
      double e;
      double fc;
      boolean_T exitg1;
      fc = fb;
      c = b;
      e = 0.0;
      d = 0.0;
      exitg1 = false;
      while ((!exitg1) && ((fb != 0.0) && (a != b))) {
        double m;
        double toler;
        if ((fb > 0.0) == (fc > 0.0)) {
          c = a;
          fc = fx;
          d = b - a;
          e = d;
        }
        if (std::abs(fc) < std::abs(fb)) {
          a = b;
          b = c;
          c = a;
          fx = fb;
          fb = fc;
          fc = fx;
        }
        m = 0.5 * (c - b);
        toler = 4.4408920985006262E-16 * std::fmax(std::abs(b), 1.0);
        if ((std::abs(m) <= toler) || (fb == 0.0)) {
          exitg1 = true;
        } else {
          if ((std::abs(e) < toler) || (std::abs(fx) <= std::abs(fb))) {
            d = m;
            e = m;
          } else {
            double s;
            s = fb / fx;
            if (a == c) {
              dx = 2.0 * m * s;
              a_eff = 1.0 - s;
            } else {
              a_eff = fx / fc;
              r = fb / fc;
              dx = s * (2.0 * m * a_eff * (a_eff - r) - (b - a) * (r - 1.0));
              a_eff = (a_eff - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (dx > 0.0) {
              a_eff = -a_eff;
            } else {
              dx = -dx;
            }
            if ((2.0 * dx < 3.0 * m * a_eff - std::abs(toler * a_eff)) &&
                (dx < std::abs(0.5 * e * a_eff))) {
              e = d;
              d = dx / a_eff;
            } else {
              d = m;
              e = m;
            }
          }
          a = b;
          fx = fb;
          if (std::abs(d) > toler) {
            b += d;
          } else if (b > c) {
            b -= toler;
          } else {
            b += toler;
          }
          a_eff = b - varargin_3;
          r = std::sin(a_eff);
          a_eff = std::cos(a_eff);
          fb = (2.59283849675502 * r * (a_eff * a_eff) +
                3.1415926535897931 * std::abs(r) * r * a_eff) -
               varargin_4;
        }
      }
    }
  }
  return b;
}

//
// Arguments    : double varargin_3
//                double varargin_4
// Return Type  : double
//
double g_fzero(double varargin_3, double varargin_4)
{
  double a_eff;
  double b;
  double fx;
  double r;
  r = std::sin(0.58423967555431222 - varargin_3);
  a_eff = std::cos(0.58423967555431222 - varargin_3);
  fx = (0.35596863975551224 * r * (a_eff * a_eff) +
        3.1415926535897931 * std::abs(r) * r * a_eff) -
       varargin_4;
  if (fx == 0.0) {
    b = 0.58423967555431222;
  } else {
    double a;
    double dx;
    double fb;
    int exitg2;
    dx = 0.011684793511086245;
    a = 0.58423967555431222;
    b = 0.58423967555431222;
    fb = fx;
    do {
      exitg2 = 0;
      if ((fx > 0.0) == (fb > 0.0)) {
        dx *= 1.4142135623730951;
        a = 0.58423967555431222 - dx;
        a_eff = (0.58423967555431222 - dx) - varargin_3;
        r = std::sin(a_eff);
        a_eff = std::cos(a_eff);
        fx = (0.35596863975551224 * r * (a_eff * a_eff) +
              3.1415926535897931 * std::abs(r) * r * a_eff) -
             varargin_4;
        if (std::isinf(fx) || std::isnan(fx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if (std::isinf(0.58423967555431222 - dx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if ((fx > 0.0) != (fb > 0.0)) {
          exitg2 = 2;
        } else {
          b = dx + 0.58423967555431222;
          a_eff = (dx + 0.58423967555431222) - varargin_3;
          r = std::sin(a_eff);
          a_eff = std::cos(a_eff);
          fb = (0.35596863975551224 * r * (a_eff * a_eff) +
                3.1415926535897931 * std::abs(r) * r * a_eff) -
               varargin_4;
          if (std::isinf(fb) || std::isnan(fb)) {
            b = rtNaN;
            exitg2 = 1;
          } else if (std::isinf(dx + 0.58423967555431222)) {
            b = rtNaN;
            exitg2 = 1;
          }
        }
      } else {
        exitg2 = 2;
      }
    } while (exitg2 == 0);
    if (exitg2 != 1) {
      double c;
      double d;
      double e;
      double fc;
      boolean_T exitg1;
      fc = fb;
      c = b;
      e = 0.0;
      d = 0.0;
      exitg1 = false;
      while ((!exitg1) && ((fb != 0.0) && (a != b))) {
        double m;
        double toler;
        if ((fb > 0.0) == (fc > 0.0)) {
          c = a;
          fc = fx;
          d = b - a;
          e = d;
        }
        if (std::abs(fc) < std::abs(fb)) {
          a = b;
          b = c;
          c = a;
          fx = fb;
          fb = fc;
          fc = fx;
        }
        m = 0.5 * (c - b);
        toler = 4.4408920985006262E-16 * std::fmax(std::abs(b), 1.0);
        if ((std::abs(m) <= toler) || (fb == 0.0)) {
          exitg1 = true;
        } else {
          if ((std::abs(e) < toler) || (std::abs(fx) <= std::abs(fb))) {
            d = m;
            e = m;
          } else {
            double s;
            s = fb / fx;
            if (a == c) {
              dx = 2.0 * m * s;
              a_eff = 1.0 - s;
            } else {
              a_eff = fx / fc;
              r = fb / fc;
              dx = s * (2.0 * m * a_eff * (a_eff - r) - (b - a) * (r - 1.0));
              a_eff = (a_eff - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (dx > 0.0) {
              a_eff = -a_eff;
            } else {
              dx = -dx;
            }
            if ((2.0 * dx < 3.0 * m * a_eff - std::abs(toler * a_eff)) &&
                (dx < std::abs(0.5 * e * a_eff))) {
              e = d;
              d = dx / a_eff;
            } else {
              d = m;
              e = m;
            }
          }
          a = b;
          fx = fb;
          if (std::abs(d) > toler) {
            b += d;
          } else if (b > c) {
            b -= toler;
          } else {
            b += toler;
          }
          a_eff = b - varargin_3;
          r = std::sin(a_eff);
          a_eff = std::cos(a_eff);
          fb = (0.35596863975551224 * r * (a_eff * a_eff) +
                3.1415926535897931 * std::abs(r) * r * a_eff) -
               varargin_4;
        }
      }
    }
  }
  return b;
}

//
// Arguments    : double varargin_3
//                double varargin_4
// Return Type  : double
//
double h_fzero(double varargin_3, double varargin_4)
{
  double a_eff;
  double b;
  double fx;
  double r;
  r = std::sin(-0.58423967555431222 - varargin_3);
  a_eff = std::cos(-0.58423967555431222 - varargin_3);
  fx = (0.35596863975551224 * r * (a_eff * a_eff) +
        3.1415926535897931 * std::abs(r) * r * a_eff) -
       varargin_4;
  if (fx == 0.0) {
    b = -0.58423967555431222;
  } else {
    double a;
    double dx;
    double fb;
    int exitg2;
    dx = -0.011684793511086245;
    a = -0.58423967555431222;
    b = -0.58423967555431222;
    fb = fx;
    do {
      exitg2 = 0;
      if ((fx > 0.0) == (fb > 0.0)) {
        dx *= 1.4142135623730951;
        a = -0.58423967555431222 - dx;
        a_eff = (-0.58423967555431222 - dx) - varargin_3;
        r = std::sin(a_eff);
        a_eff = std::cos(a_eff);
        fx = (0.35596863975551224 * r * (a_eff * a_eff) +
              3.1415926535897931 * std::abs(r) * r * a_eff) -
             varargin_4;
        if (std::isinf(fx) || std::isnan(fx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if (std::isinf(-0.58423967555431222 - dx)) {
          b = rtNaN;
          exitg2 = 1;
        } else if ((fx > 0.0) != (fb > 0.0)) {
          exitg2 = 2;
        } else {
          b = dx + -0.58423967555431222;
          a_eff = (dx + -0.58423967555431222) - varargin_3;
          r = std::sin(a_eff);
          a_eff = std::cos(a_eff);
          fb = (0.35596863975551224 * r * (a_eff * a_eff) +
                3.1415926535897931 * std::abs(r) * r * a_eff) -
               varargin_4;
          if (std::isinf(fb) || std::isnan(fb)) {
            b = rtNaN;
            exitg2 = 1;
          } else if (std::isinf(dx + -0.58423967555431222)) {
            b = rtNaN;
            exitg2 = 1;
          }
        }
      } else {
        exitg2 = 2;
      }
    } while (exitg2 == 0);
    if (exitg2 != 1) {
      double c;
      double d;
      double e;
      double fc;
      boolean_T exitg1;
      fc = fb;
      c = b;
      e = 0.0;
      d = 0.0;
      exitg1 = false;
      while ((!exitg1) && ((fb != 0.0) && (a != b))) {
        double m;
        double toler;
        if ((fb > 0.0) == (fc > 0.0)) {
          c = a;
          fc = fx;
          d = b - a;
          e = d;
        }
        if (std::abs(fc) < std::abs(fb)) {
          a = b;
          b = c;
          c = a;
          fx = fb;
          fb = fc;
          fc = fx;
        }
        m = 0.5 * (c - b);
        toler = 4.4408920985006262E-16 * std::fmax(std::abs(b), 1.0);
        if ((std::abs(m) <= toler) || (fb == 0.0)) {
          exitg1 = true;
        } else {
          if ((std::abs(e) < toler) || (std::abs(fx) <= std::abs(fb))) {
            d = m;
            e = m;
          } else {
            double s;
            s = fb / fx;
            if (a == c) {
              dx = 2.0 * m * s;
              a_eff = 1.0 - s;
            } else {
              a_eff = fx / fc;
              r = fb / fc;
              dx = s * (2.0 * m * a_eff * (a_eff - r) - (b - a) * (r - 1.0));
              a_eff = (a_eff - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (dx > 0.0) {
              a_eff = -a_eff;
            } else {
              dx = -dx;
            }
            if ((2.0 * dx < 3.0 * m * a_eff - std::abs(toler * a_eff)) &&
                (dx < std::abs(0.5 * e * a_eff))) {
              e = d;
              d = dx / a_eff;
            } else {
              d = m;
              e = m;
            }
          }
          a = b;
          fx = fb;
          if (std::abs(d) > toler) {
            b += d;
          } else if (b > c) {
            b -= toler;
          } else {
            b += toler;
          }
          a_eff = b - varargin_3;
          r = std::sin(a_eff);
          a_eff = std::cos(a_eff);
          fb = (0.35596863975551224 * r * (a_eff * a_eff) +
                3.1415926535897931 * std::abs(r) * r * a_eff) -
               varargin_4;
        }
      }
    }
  }
  return b;
}

} // namespace coder

//
// File trailer for fzero.cpp
//
// [EOF]
//

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

//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: McFoamy_FM_v3_data.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 05-Aug-2021 13:41:18
//

// Include Files
#include "McFoamy_FM_v3_data.h"
#include "rt_nonfinite.h"

//
// File trailer for McFoamy_FM_v3_data.cpp
//
// [EOF]
//

//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: McFoamy_FM_v3_initialize.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 05-Aug-2021 13:41:18
//

// Include Files
#include "McFoamy_FM_v3_initialize.h"
#include "rt_nonfinite.h"

// Function Definitions
//
// Arguments    : void
// Return Type  : void
//
void McFoamy_FM_v3_initialize()
{
}

//
// File trailer for McFoamy_FM_v3_initialize.cpp
//
// [EOF]
//

//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: McFoamy_FM_v3_terminate.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 05-Aug-2021 13:41:18
//

// Include Files
#include "McFoamy_FM_v3_terminate.h"
#include "rt_nonfinite.h"

// Function Definitions
//
// Arguments    : void
// Return Type  : void
//
void McFoamy_FM_v3_terminate()
{
  // (no terminate code required)
}

//
// File trailer for McFoamy_FM_v3_terminate.cpp
//
// [EOF]
//

//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: rt_nonfinite.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 05-Aug-2021 13:41:18
//

// Abstract:
//      MATLAB for code generation function to initialize non-finites,
//      (Inf, NaN and -Inf).
// Include Files
#include "rt_nonfinite.h"
#include <cmath>
#include <limits>

real_T rtNaN{std::numeric_limits<real_T>::quiet_NaN()};
real_T rtInf{std::numeric_limits<real_T>::infinity()};
real_T rtMinusInf{-std::numeric_limits<real_T>::infinity()};
real32_T rtNaNF{std::numeric_limits<real32_T>::quiet_NaN()};
real32_T rtInfF{std::numeric_limits<real32_T>::infinity()};
real32_T rtMinusInfF{-std::numeric_limits<real32_T>::infinity()};

//
// File trailer for rt_nonfinite.cpp
//
// [EOF]
//
//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: rtGetInf.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 05-Aug-2021 13:41:18
//

// Abstract:
//       MATLAB for code generation function to initialize non-finite, Inf and
//       MinusInf
// Include Files
#include "rtGetInf.h"
#include "rt_nonfinite.h"

// Function: rtGetInf
// ==================================================================
//  Abstract:
//  Initialize rtInf needed by the generated code.
real_T rtGetInf(void)
{
  return rtInf;
}

// Function: rtGetInfF
// =================================================================
//  Abstract:
//  Initialize rtInfF needed by the generated code.
real32_T rtGetInfF(void)
{
  return rtInfF;
}

// Function: rtGetMinusInf
// =============================================================
//  Abstract:
//  Initialize rtMinusInf needed by the generated code.
real_T rtGetMinusInf(void)
{
  return rtMinusInf;
}

// Function: rtGetMinusInfF
// ============================================================
//  Abstract:
//  Initialize rtMinusInfF needed by the generated code.
real32_T rtGetMinusInfF(void)
{
  return rtMinusInfF;
}

//
// File trailer for rtGetInf.cpp
//
// [EOF]
//
//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: rtGetNaN.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 05-Aug-2021 13:41:18
//

// Abstract:
//       MATLAB for code generation function to initialize non-finite, NaN
// Include Files
#include "rtGetNaN.h"
#include "rt_nonfinite.h"

// Function: rtGetNaN
// ======================================================================
//  Abstract:
// Initialize rtNaN needed by the generated code.
//  NaN is initialized as non-signaling. Assumes IEEE.
real_T rtGetNaN(void)
{
  return rtNaN;
}

// Function: rtGetNaNF
// =====================================================================
//  Abstract:
//  Initialize rtNaNF needed by the generated code.
//  NaN is initialized as non-signaling. Assumes IEEE
real32_T rtGetNaNF(void)
{
  return rtNaNF;
}

//
// File trailer for rtGetNaN.cpp
//
// [EOF]
//
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
