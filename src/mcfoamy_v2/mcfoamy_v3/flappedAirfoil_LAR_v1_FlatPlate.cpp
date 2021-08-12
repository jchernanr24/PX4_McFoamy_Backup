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
