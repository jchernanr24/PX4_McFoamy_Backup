//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: McFoamy_FM_v2.cpp
//
// MATLAB Coder version            : 4.2
// C/C++ source code generated on  : 09-Aug-2019 16:50:12
//

// Include Files

#include <string.h>
#include <math.h>

//#include "mcfoamy_v2/McFoamy_FM_v2.h"
// Function Declarations
static void McFoamy_Airfoil_Simplified(double alp, double cf, double c, double
  def, double AR, double Cd0, double Cd90, double alp0, double *CN, double *CL,
  double *CD, double *CM);
static void b_asin(double *x);
static void b_atan2(const double y[7], const double x[7], double r[7]);
static int b_bsearch(const double x[10], double xi);
static double b_interp1(const double varargin_1[10], const double varargin_2[10],
  double varargin_3);
static double b_norm(const double x[3]);
static void b_power(const double a[3], double y[3]);
static void b_sign(double x[4]);
static void b_sqrt(double *x);
static void c_atan2(const double y[3], const double x[3], double r[3]);
static double c_interp1(const double varargin_1[2], const double varargin_2[2],
  double varargin_3);
static void c_power(const double a[4], double y[4]);
static void c_sign(double *x);
static void c_sqrt(double x[7]);
static void cross(const double a[3], const double b[3], double c[3]);
static void d_atan2(const double y[4], const double x[4], double r[4]);
static void d_sqrt(double x[3]);
static void e_sqrt(double x[4]);
static double interp1(const double varargin_1[10], const double varargin_2[10],
                      double varargin_3);
static double interp2_dispatch(double V[10][11], double Xq, double Yq, double
  extrapval);
static void power(const double a[7], double y[7]);
static double rt_atan2d_snf(double u0, double u1);
static double rt_powd_snf(double u0, double u1);

// Function Definitions

//
// Arguments    : double alp
//                double cf
//                double c
//                double def
//                double AR
//                double Cd0
//                double Cd90
//                double alp0
//                double *CN
//                double *CL
//                double *CD
//                double *CM
// Return Type  : void
//
static void McFoamy_Airfoil_Simplified(double alp, double cf, double c, double
  def, double AR, double Cd0, double Cd90, double alp0, double *CN, double *CL,
  double *CD, double *CM)
{
  double LowAlpEnd;
  static const double dv20[10] = { 0.1666, 0.333, 0.4, 0.5, 1.0, 1.25, 2.0, 3.0,
    4.0, 6.0 };

  static const double dv21[10] = { 0.55850536063818546, 0.62831853071795862,
    0.59341194567807209, 0.66322511575784515, 0.55850536063818546,
    0.36651914291880922, 0.27925268031909273, 0.19198621771937624,
    0.17453292519943295, 0.13962634015954636 };

  double HighAlpStart;
  static const double dv22[10] = { 0.69813170079773179, 1.0471975511965976,
    0.95993108859688125, 0.97738438111682457, 0.69813170079773179,
    0.50614548307835561, 0.48869219055841229, 0.41887902047863906,
    0.38397243543875248, 0.3490658503988659 };

  double CLAlp;
  double d14;
  double alp0_eff;
  double LowAlpEnd_N;
  double LowAlpStart_P;
  double LowAlpStart_N;
  double HighAlpEnd_P;
  double HighAlpStart_N;
  double HighAlpEnd_N;
  double b_gamma;
  double CT;
  double a_eff;
  double x;
  double a;
  double x_tmp;
  double y;
  double b_x;
  double c_x;
  double b_y;
  double d_x;
  double CN1;
  double b_CN1[2];
  double b_LowAlpStart_P[2];
  double b_CLAlp[2];

  //  Slightly different from the original McFoamy_Airfoil_Simplified
  //  Should only be used for MATLAB Code to get MEX file.
  // ----------------- STALL ANGLES AND HIGH AOA ANGLES -----------------%
  //  AR_pre = [0.5,0.75,1,1.25,1.5,1.75,2, 3,4,6];
  //  alpStall_pre = [44,37,34,22,20,18,18, 18,17,17]*pi/180;
  //  alpLim_pre = [35,33,28,20,15,14,13, 13,13,13]*pi/180; % By Torres
  // alpLim_pre = [32,36,34,38,32,21,16,11,10,8]*pi/180; % Alp_Limit
  //  Alp_Limit
  LowAlpEnd = b_interp1(dv20, dv21, AR);
  HighAlpStart = b_interp1(dv20, dv22, AR);

  // ----------------------- LIFT CURVE SLOPE ---------------------------%
  //  From McCormick
  CLAlp = 6.2831853071795862 * AR / (AR + 2.0 * (AR + 4.0) / (AR + 2.0));

  //  -------- POTENTIAL LIFT AND VORTEX LIFT COEFFICIENTS --------------%
  // --------------------- FLAP DEFLECTION ------------------------------%
  d14 = cf / c;
  if (d14 == 1.0) {
    alp0_eff = alp0;
    LowAlpEnd_N = -LowAlpEnd;
    LowAlpStart_P = -LowAlpEnd + M_PI;
    LowAlpStart_N = LowAlpEnd - M_PI;
    HighAlpEnd_P = M_PI - HighAlpStart;
    HighAlpStart_N = -HighAlpStart;
    HighAlpEnd_N = -(M_PI - HighAlpStart);
    alp += def;
    b_gamma = 0.0;
  } else {
    // eta_pre = [0.81*0.7 0.8*0.7 0.685*0.7 0.535 0.455 0.405 0.37 0.345];
    //  For Linear CL vs Alp curve
    //  alp0_eff1 = alp0 - delCL/CLAlp;
    //  CLmaxP1 = CLAlp*(alpStallP - alp0) + delCLmax;
    //  CLmaxN1 = CLAlp*(alpStallN - alp0) + delCLmax;
    //  alpStallP_eff1 = alp0_eff1 + CLmaxP1/CLAlp;
    //  alpStallN_eff1 = alp0_eff1 + CLmaxN1/CLAlp;
    //  For Non-Linear CL vs Alp curve
    // fun = @(x) Kp*sin(0-x)*(cos(0-x))^2 + Kv*abs(sin(0-x))*sin(0-x)*cos(0-x) - delCL;
    //      alp0_eff = fzero(@myfun_alp0_eff,x0,[],Kp,Kv,delCL);
    // fun = @(x) Kp*sin(x-alp0_eff)*(cos(x-alp0_eff))^2 + Kv*abs(sin(x-alp0_eff))*sin(x-alp0_eff)*cos(x-alp0_eff) - CLmaxP;
    //      LowAlpEnd_P = fzero(@myfun_alpStall_eff,x0,[],Kp,Kv,alp0_eff,CLmaxP);
    // alpStallP_eff = fzero(@myfun_alpStall_eff,x0,[],Kp,Kv,alp0_eff,CLmaxP);
    // fun = @(x) Kp*sin(x-alp0_eff)*(cos(x-alp0_eff))^2 + Kv*abs(sin(x-alp0_eff))*sin(x-alp0_eff)*cos(x-alp0_eff) - CLmaxN;
    //      LowAlpEnd_N = fzero(@myfun_alpStall_eff,x0,[],Kp,Kv,alp0_eff,CLmaxN);
    // alpStallN_eff = fzero(@myfun_alpStall_eff,x0,[],Kp,Kv,alp0_eff,CLmaxN);
    if (fabs(AR - 2.092) < 0.0009) {
      // wing
      CT = rt_powd_snf(def, 3.0);
      HighAlpEnd_P = def * def;
      alp0_eff = (0.1145 * CT + 0.0 * HighAlpEnd_P) + -0.3372 * def;
      CT *= 0.0678;
      LowAlpEnd = ((CT + 0.0009 * HighAlpEnd_P) + -0.1636 * def) + 0.2713;
      LowAlpEnd_N = ((CT + -0.0009 * HighAlpEnd_P) + -0.1636 * def) + -0.2713;
    } else if (fabs(AR - 0.227938) < 0.0009) {
      // body
      alp0_eff = (0.0 * rt_powd_snf(def, 3.0) + 0.0 * (def * def)) + 0.0 * def;
      LowAlpEnd = alp0_eff + 0.5842;
      LowAlpEnd_N = alp0_eff + -0.5842;
    } else if (fabs(AR - 1.298412) < 0.0009) {
      // rud
      alp0_eff = (0.1263 * rt_powd_snf(def, 3.0) + 0.0 * (def * def)) + -0.3456 *
        def;
      CT = 0.0951 * rt_powd_snf(def, 3.0);
      LowAlpEnd = ((CT + 0.0015 * (def * def)) + -0.2331 * def) + 0.361;
      LowAlpEnd_N = ((CT + -0.0015 * (def * def)) + -0.2331 * def) + -0.361;
    } else {
      // tail
      alp0_eff = (0.1381 * rt_powd_snf(def, 3.0) + -0.0 * (def * def)) + -0.3734
        * def;
      CT = 0.1196 * rt_powd_snf(def, 3.0);
      LowAlpEnd = ((CT + 0.0006 * (def * def)) + -0.3079 * def) + 0.3726;
      LowAlpEnd_N = ((CT + -0.0006 * (def * def)) + -0.3079 * def) + -0.3726;
    }

    LowAlpStart_P = LowAlpEnd_N + M_PI;
    LowAlpStart_N = LowAlpEnd - M_PI;

    //  High AoA parameters using equivalent flat plate method
    CT = c - cf;
    HighAlpEnd_P = (CT * CT + cf * cf) + 2.0 * CT * cf * cosf(fabs(def));
    b_sqrt(&HighAlpEnd_P);
    b_gamma = sinf(def) * cf / HighAlpEnd_P;
    b_asin(&b_gamma);
    HighAlpEnd_P = M_PI - HighAlpStart;
    HighAlpStart_N = -HighAlpStart;
    HighAlpEnd_N = -(M_PI - HighAlpStart);
  }

  //     %% Calculate Aerodynamic coefficients based on regime
  if ((alp >= -M_PI) && (alp <= LowAlpStart_N)) {
    LowAlpStart_P = alp - alp0_eff;
    CT = cosf(LowAlpStart_P + M_PI);
    alp0_eff = sinf(LowAlpStart_P + M_PI);
    *CL = CLAlp * alp0_eff * (CT * CT) + M_PI * fabs(alp0_eff)
      * alp0_eff * CT;
    CT = sinf(LowAlpStart_P + M_PI);
    HighAlpEnd_P = cosf(LowAlpStart_P + M_PI);
    alp0_eff = fabs(CT);
    *CD = Cd0 + fabs(CLAlp * alp0_eff * CT * HighAlpEnd_P +
                         M_PI * rt_powd_snf(CT, 3.0));
    *CN = -(CLAlp * CT * fabs(HighAlpEnd_P) + M_PI * alp0_eff *
            CT);

    // CM = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
    *CM = -(0.25 - 0.175 * (1.0 - fabs((LowAlpStart_P + M_PI)
              - M_PI) / M_PI_2)) * *CN;

    //  ------------------------------------------------- %
  } else if ((alp > LowAlpStart_N) && (alp < HighAlpEnd_N)) {
    //  Linear End
    a_eff = (LowAlpStart_N - alp0_eff) + M_PI;
    x = sinf(a_eff);
    a = cosf(a_eff);
    y = fabs(sinf(a_eff));
    b_x = sinf(a_eff);
    x_tmp = cosf(a_eff);
    HighAlpStart_N = fabs(sinf(a_eff));
    c_x = sinf(a_eff);
    d_x = cosf(a_eff);
    alp0_eff = rt_powd_snf(sinf(a_eff), 3.0);
    CN1 = -(CLAlp * sinf(a_eff) * cosf(a_eff) + M_PI * fabs(sinf(a_eff)) * sinf(a_eff));

    // CM1 = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
    // CD1 = Cd0 + abs(CL1)*abs(tan(a_eff));
    // CN1 = -sqrt(CL1^2 + CD1^2);
    LowAlpStart_P = fabs(a_eff - M_PI);

    //  NonLinear End
    a_eff = HighAlpEnd_N + b_gamma;

    // - alp0_eff;
    //  Determining effective Cd90
    if ((d14 == 1.0) || (d14 == 0.0)) {
      b_gamma = Cd90;
    } else {
      d14 = HighAlpEnd_N;
      c_sign(&d14);
      if (d14 == 1.0) {
        CT = def * 180.0 / M_PI;
        b_gamma = (-1.2963E-5 * (CT * CT) + 0.00361111 * (def * 180.0 /
                    M_PI)) + Cd90;
      } else {
        CT = def * 180.0 / M_PI;
        b_gamma = (-1.2963E-5 * (CT * CT) - 0.00361111 * (def * 180.0 /
                    M_PI)) + Cd90;
      }
    }

    //  Hoerner Model with Lindenberg Correction
    CT = sinf(fabs(a_eff));
    b_gamma = -b_gamma * (1.0 / (0.56 + 0.44 * CT) - 0.41 * (1.0 - expf
      (-17.0 / AR))) * CT;
    HighAlpEnd_P = 0.5 * Cd0 * cosf(a_eff);

    //  Interpolation
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    b_CN1[0] = CN1;
    b_CN1[1] = b_gamma;
    *CN = c_interp1(b_LowAlpStart_P, b_CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    b_CLAlp[0] = CLAlp * x * (a * a) + M_PI * y * b_x * x_tmp;
    CT = sinf(a_eff);
    b_CLAlp[1] = b_gamma * cosf(a_eff) - HighAlpEnd_P * CT;
    *CL = c_interp1(b_LowAlpStart_P, b_CLAlp, alp);
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    b_CLAlp[0] = Cd0 + fabs(CLAlp * HighAlpStart_N * c_x * d_x +
      M_PI * alp0_eff);
    b_CLAlp[1] = b_gamma * CT + HighAlpEnd_P * cosf(a_eff);
    *CD = c_interp1(b_LowAlpStart_P, b_CLAlp, alp);
    b_LowAlpStart_P[0] = LowAlpStart_N;
    b_LowAlpStart_P[1] = HighAlpEnd_N;
    b_CLAlp[0] = -(0.25 - 0.175 * (1.0 - LowAlpStart_P / M_PI_2)) *
      CN1;
    b_CLAlp[1] = -(0.25 - 0.175 * (1.0 - fabs(a_eff) / M_PI_2)) *
      b_gamma;
    *CM = c_interp1(b_LowAlpStart_P, b_CLAlp, alp);

    //  ------------------------------------------------- %
  } else if ((alp >= HighAlpEnd_N) && (alp <= HighAlpStart_N)) {
    //  Equivalent flat plate
    LowAlpStart_P = alp + b_gamma;

    // - alp0_eff;
    //  Determining effective Cd90
    if ((d14 == 1.0) || (d14 == 0.0)) {
      b_gamma = Cd90;
    } else {
      d14 = alp;
      c_sign(&d14);
      if (d14 == 1.0) {
        a = def * 180.0 / M_PI;
        b_gamma = (-1.2963E-5 * (a * a) + 0.00361111 * (def * 180.0 /
                    M_PI)) + Cd90;
      } else {
        CT = def * 180.0 / M_PI;
        b_gamma = (-1.2963E-5 * (CT * CT) - 0.00361111 * CT) + Cd90;
      }
    }

    //  Hoerner Model
    CT = sinf(fabs(LowAlpStart_P));
    *CN = -b_gamma * (1.0 / (0.56 + 0.44 * CT) - 0.41 * (1.0 - expf(-17.0 /
      AR))) * CT;
    CT = 0.5 * Cd0 * cosf(LowAlpStart_P);
    *CL = *CN * cosf(LowAlpStart_P) - CT * sinf(LowAlpStart_P);
    *CD = *CN * sinf(LowAlpStart_P) + CT * cosf(LowAlpStart_P);
    *CM = -(0.25 - 0.175 * (1.0 - fabs(LowAlpStart_P) / M_PI_2))
      * *CN;

    //  ------------------------------------------------- %
  } else if ((alp > HighAlpStart_N) && (alp < LowAlpEnd_N)) {
    //  Linear End
    a_eff = LowAlpEnd_N - alp0_eff;
    x_tmp = sinf(a_eff);
    a = cosf(a_eff);
    y = fabs(sinf(a_eff));
    x = sinf(a_eff);
    b_x = cosf(a_eff);
    b_y = fabs(sinf(a_eff));
    c_x = sinf(a_eff);
    LowAlpStart_P = cosf(a_eff);
    CN1 = rt_powd_snf(sinf(a_eff), 3.0);
    d_x = sinf(a_eff);
    HighAlpEnd_N = fabs(sinf(a_eff));
    LowAlpStart_N = sinf(a_eff);
    alp0_eff = fabs(sinf(a_eff));
    HighAlpEnd_P = sinf(a_eff);

    //  NonLinear End
    a_eff = HighAlpStart_N + b_gamma;

    // - alp0_eff;
    if ((d14 == 1.0) || (d14 == 0.0)) {
      b_gamma = Cd90;
    } else {
      d14 = HighAlpStart_N;
      c_sign(&d14);
      if (d14 == 1.0) {
        CT = def * 180.0 / M_PI;
        b_gamma = (-1.2963E-5 * (CT * CT) + 0.00361111 * (def * 180.0 /
                    M_PI)) + Cd90;
      } else {
        CT = def * 180.0 / M_PI;
        b_gamma = (-1.2963E-5 * (CT * CT) - 0.00361111 * (def * 180.0 /
                    M_PI)) + Cd90;
      }
    }

    CT = sinf(fabs(a_eff));
    b_gamma = -b_gamma * (1.0 / (0.56 + 0.44 * CT) - 0.41 * (1.0 - expf
      (-17.0 / AR))) * CT;
    CT = 0.5 * Cd0 * cosf(a_eff);

    //  Interpolation
    b_CN1[0] = LowAlpEnd_N;
    b_CN1[1] = HighAlpStart_N;
    b_CLAlp[0] = CLAlp * d_x * LowAlpStart_P + M_PI * HighAlpEnd_N
      * LowAlpStart_N;
    b_CLAlp[1] = b_gamma;
    *CN = c_interp1(b_CN1, b_CLAlp, alp);
    b_CN1[0] = LowAlpEnd_N;
    b_CN1[1] = HighAlpStart_N;
    b_CLAlp[0] = CLAlp * x_tmp * (a * a) + M_PI * y * x * b_x;
    b_CLAlp[1] = b_gamma * cosf(a_eff) - CT * sinf(a_eff);
    *CL = c_interp1(b_CN1, b_CLAlp, alp);
    b_CN1[0] = LowAlpEnd_N;
    b_CN1[1] = HighAlpStart_N;
    b_CLAlp[0] = Cd0 + fabs(CLAlp * b_y * c_x * LowAlpStart_P +
      M_PI * CN1);
    b_CLAlp[1] = b_gamma * sinf(a_eff) + CT * cosf(a_eff);
    *CD = c_interp1(b_CN1, b_CLAlp, alp);
    b_CN1[0] = LowAlpEnd_N;
    b_CN1[1] = HighAlpStart_N;
    b_CLAlp[0] = -0.53407075111026481 * alp0_eff * HighAlpEnd_P;
    b_CLAlp[1] = -(0.25 - 0.175 * (1.0 - fabs(a_eff) / M_PI_2)) *
      b_gamma;
    *CM = c_interp1(b_CN1, b_CLAlp, alp);

    //  ------------------------------------------------- %
  } else if ((alp >= LowAlpEnd_N) && (alp <= LowAlpEnd)) {
    d14 = alp - alp0_eff;
    a = cosf(d14);
    *CL = CLAlp * sinf(d14) * (a * a) + M_PI * fabs(sinf(d14)) * sinf(d14) * cosf(d14);
    CT = sinf(d14);
    HighAlpEnd_P = cosf(d14);
    alp0_eff = fabs(CT);
    *CD = (Cd0 + fabs(CLAlp * alp0_eff * CT * HighAlpEnd_P)) +
      M_PI * fabs(rt_powd_snf(CT, 3.0));
    *CN = CLAlp * CT * HighAlpEnd_P + M_PI * alp0_eff * CT;
    CT = sinf(d14);
    *CM = -0.53407075111026481 * fabs(CT) * CT;

    //  ------------------------------------------------- %
  } else if ((alp > LowAlpEnd) && (alp < HighAlpStart)) {
    //  Linear End
    a_eff = LowAlpEnd - alp0_eff;
    x = sinf(a_eff);
    CT = cosf(a_eff);
    HighAlpStart_N = sinf(a_eff);
    x_tmp = sinf(a_eff);
    b_x = cosf(a_eff);
    y = fabs(sinf(a_eff));
    c_x = cosf(a_eff);
    CN1 = sinf(a_eff);
    d_x = cosf(a_eff);
    b_y = fabs(sinf(a_eff));
    HighAlpEnd_N = fabs(sinf(a_eff));
    LowAlpStart_N = sinf(a_eff);

    //  NonLinear End
    a_eff = HighAlpStart + b_gamma;

    // - alp0_eff;
    if ((d14 == 1.0) || (d14 == 0.0)) {
      b_gamma = Cd90;
    } else {
      d14 = HighAlpStart;
      c_sign(&d14);
      if (d14 == 1.0) {
        a = def * 180.0 / M_PI;
        b_gamma = (-1.2963E-5 * (a * a) + 0.00361111 * (def * 180.0 /
                    M_PI)) + Cd90;
      } else {
        a = def * 180.0 / M_PI;
        b_gamma = (-1.2963E-5 * (a * a) - 0.00361111 * (def * 180.0 /
                    M_PI)) + Cd90;
      }
    }

    b_gamma = b_gamma * (1.0 / (0.56 + 0.44 * sinf(fabs(a_eff))) - 0.41 *
                         (1.0 - expf(-17.0 / AR))) * sinf(fabs(a_eff));
    d14 = 0.5 * Cd0 * cosf(a_eff);

    //  Interpolation
    b_CN1[0] = LowAlpEnd;
    b_CN1[1] = HighAlpStart;
    b_CLAlp[0] = CLAlp * CN1 * d_x + M_PI * b_y * CN1;
    b_CLAlp[1] = b_gamma;
    *CN = c_interp1(b_CN1, b_CLAlp, alp);
    b_CN1[0] = LowAlpEnd;
    b_CN1[1] = HighAlpStart;
    b_CLAlp[0] = CLAlp * x * (CT * CT) + M_PI * fabs
      (HighAlpStart_N) * x_tmp * b_x;
    b_CLAlp[1] = b_gamma * cosf(a_eff) - d14 * sinf(a_eff);
    *CL = c_interp1(b_CN1, b_CLAlp, alp);
    b_CN1[0] = LowAlpEnd;
    b_CN1[1] = HighAlpStart;
    b_CLAlp[0] = Cd0 + fabs(CLAlp * y * x_tmp * c_x + M_PI *
      rt_powd_snf(CN1, 3.0));
    b_CLAlp[1] = b_gamma * sinf(a_eff) + d14 * cosf(a_eff);
    *CD = c_interp1(b_CN1, b_CLAlp, alp);
    b_CN1[0] = LowAlpEnd;
    b_CN1[1] = HighAlpStart;
    b_CLAlp[0] = -0.53407075111026481 * HighAlpEnd_N * LowAlpStart_N;
    b_CLAlp[1] = -(0.25 - 0.175 * (1.0 - fabs(a_eff) / M_PI_2)) *
      b_gamma;
    *CM = c_interp1(b_CN1, b_CLAlp, alp);

    //  ------------------------------------------------- %
  } else if ((alp >= HighAlpStart) && (alp <= HighAlpEnd_P)) {
    //  Equivalent flat plate
    HighAlpEnd_P = alp + b_gamma;

    // - alp0_eff;
    //  Determining effective Cd90
    if ((d14 == 1.0) || (d14 == 0.0)) {
      b_gamma = Cd90;
    } else {
      d14 = alp;
      c_sign(&d14);
      if (d14 == 1.0) {
        a = def * 180.0 / M_PI;
        b_gamma = (-1.2963E-5 * (a * a) + 0.00361111 * (def * 180.0 /
                    M_PI)) + Cd90;
      } else {
        a = def * 180.0 / M_PI;
        b_gamma = (-1.2963E-5 * (a * a) - 0.00361111 * (def * 180.0 /
                    M_PI)) + Cd90;
      }
    }

    //  Hoerner Model
    CT = sinf(fabs(HighAlpEnd_P));
    *CN = b_gamma * (1.0 / (0.56 + 0.44 * CT) - 0.41 * (1.0 - expf(-17.0 /
      AR))) * CT;
    CT = 0.5 * Cd0 * cosf(HighAlpEnd_P);
    *CL = *CN * cosf(HighAlpEnd_P) - CT * sinf(HighAlpEnd_P);
    *CD = *CN * sinf(HighAlpEnd_P) + CT * cosf(HighAlpEnd_P);
    *CM = -(0.25 - 0.175 * (1.0 - fabs(HighAlpEnd_P) / M_PI_2)) *
      *CN;

    //  ------------------------------------------------- %
  } else if ((alp > HighAlpEnd_P) && (alp < LowAlpStart_P)) {
    //  Linear End
    a_eff = (LowAlpStart_P - alp0_eff) - M_PI;
    x = sinf(a_eff);
    a = cosf(a_eff);
    y = fabs(sinf(a_eff));
    b_x = sinf(a_eff);
    x_tmp = cosf(a_eff);
    HighAlpStart_N = fabs(sinf(a_eff));
    c_x = sinf(a_eff);
    alp0_eff = rt_powd_snf(sinf(a_eff), 3.0);
    CT = sinf(a_eff);
    CN1 = -(CLAlp * CT * cosf(a_eff) + M_PI * fabs(CT) *
            CT);
    b_y = fabs(a_eff - M_PI);

    // CM1 = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
    //  NonLinear End
    a_eff = HighAlpEnd_P + b_gamma;

    // - alp0_eff;
    //  Determining effective Cd90
    if ((d14 == 1.0) || (d14 == 0.0)) {
      b_gamma = Cd90;
    } else {
      d14 = HighAlpEnd_P;
      c_sign(&d14);
      if (d14 == 1.0) {
        CT = def * 180.0 / M_PI;
        b_gamma = (-1.2963E-5 * (CT * CT) + 0.00361111 * (def * 180.0 /
                    M_PI)) + Cd90;
      } else {
        CT = def * 180.0 / M_PI;
        b_gamma = (-1.2963E-5 * (CT * CT) - 0.00361111 * (def * 180.0 /
                    M_PI)) + Cd90;
      }
    }

    //  Hoerner Model
    b_gamma = b_gamma * (1.0 / (0.56 + 0.44 * sinf(fabs(a_eff))) - 0.41 *
                         (1.0 - expf(-17.0 / AR))) * sinf(fabs(a_eff));
    d14 = 0.5 * Cd0 * cosf(a_eff);

    //  Interpolation
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = HighAlpEnd_P;
    b_CN1[0] = CN1;
    b_CN1[1] = b_gamma;
    *CN = c_interp1(b_LowAlpStart_P, b_CN1, alp);
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = HighAlpEnd_P;
    b_CLAlp[0] = CLAlp * x * (a * a) + M_PI * y * b_x * x_tmp;
    b_CLAlp[1] = b_gamma * cosf(a_eff) - d14 * sinf(a_eff);
    *CL = c_interp1(b_LowAlpStart_P, b_CLAlp, alp);
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = HighAlpEnd_P;
    b_CLAlp[0] = Cd0 + fabs(CLAlp * HighAlpStart_N * c_x * x_tmp +
      M_PI * alp0_eff);
    b_CLAlp[1] = b_gamma * sinf(a_eff) + d14 * cosf(a_eff);
    *CD = c_interp1(b_LowAlpStart_P, b_CLAlp, alp);
    b_LowAlpStart_P[0] = LowAlpStart_P;
    b_LowAlpStart_P[1] = HighAlpEnd_P;
    b_CLAlp[0] = -(0.25 - 0.175 * (1.0 - b_y / M_PI_2)) * CN1;
    b_CLAlp[1] = -(0.25 - 0.175 * (1.0 - fabs(a_eff) / M_PI_2)) *
      b_gamma;
    *CM = c_interp1(b_LowAlpStart_P, b_CLAlp, alp);

    //  ------------------------------------------------- %
  } else {
    d14 = alp - alp0_eff;
    a_eff = (alp - alp0_eff) - M_PI;
    a = cosf(d14 - M_PI);
    alp0_eff = fabs(sinf(a_eff));
    CT = CLAlp * sinf(a_eff);
    HighAlpEnd_P = M_PI * alp0_eff * sinf(a_eff);
    *CL = CT * (a * a) + HighAlpEnd_P * cosf(a_eff);
    *CD = Cd0 + fabs(CLAlp * alp0_eff * sinf(a_eff) * cosf(a_eff) +
                         M_PI * rt_powd_snf(sinf(a_eff), 3.0));
    *CN = -(CT * cosf(a_eff) + HighAlpEnd_P);
    *CM = -(0.25 - 0.175 * (1.0 - fabs((d14 - M_PI) -
              M_PI) / M_PI_2)) * *CN;

    // CM = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
  }
}

//
// Arguments    : double *x
// Return Type  : void
//
static void b_asin(double *x)
{
  *x = asinf(*x);
}

//
// Arguments    : const double y[7]
//                const double x[7]
//                double r[7]
// Return Type  : void
//
static void b_atan2(const double y[7], const double x[7], double r[7])
{
  int k;
  for (k = 0; k < 7; k++) {
    r[k] = rt_atan2d_snf(y[k], x[k]);
  }
}

//
// Arguments    : const double x[10]
//                double xi
// Return Type  : int
//
static int b_bsearch(const double x[10], double xi)
{
  int n;
  int low_ip1;
  int high_i;
  int mid_i;
  n = 1;
  low_ip1 = 2;
  high_i = 10;
  while (high_i > low_ip1) {
    mid_i = (n + high_i) >> 1;
    if (xi >= x[mid_i - 1]) {
      n = mid_i;
      low_ip1 = mid_i + 1;
    } else {
      high_i = mid_i;
    }
  }

  return n;
}

//
// Arguments    : const double varargin_1[10]
//                const double varargin_2[10]
//                double varargin_3
// Return Type  : double
//
static double b_interp1(const double varargin_1[10], const double varargin_2[10],
  double varargin_3)
{
  double Vq;
  double y[10];
  double x[10];
  int n;
  int exitg1;
  double xtmp;
  memcpy(&y[0], &varargin_2[0], 10U * sizeof(double));
  memcpy(&x[0], &varargin_1[0], 10U * sizeof(double));
  n = 0;
  do {
    exitg1 = 0;
    if (n < 10) {
      if (isnan(varargin_1[n])) {
        exitg1 = 1;
      } else {
        n++;
      }
    } else {
      if (varargin_1[1] < varargin_1[0]) {
        for (n = 0; n < 5; n++) {
          xtmp = x[n];
          x[n] = x[9 - n];
          x[9 - n] = xtmp;
          xtmp = y[n];
          y[n] = y[9 - n];
          y[9 - n] = xtmp;
        }
      }

      Vq = NAN;
      if ((!isnan(varargin_3)) && (!(varargin_3 > x[9])) && (!(varargin_3 < x
            [0]))) {
        n = b_bsearch(x, varargin_3) - 1;
        xtmp = (varargin_3 - x[n]) / (x[n + 1] - x[n]);
        if (xtmp == 0.0) {
          Vq = y[n];
        } else if (xtmp == 1.0) {
          Vq = y[n + 1];
        } else if (y[n] == y[n + 1]) {
          Vq = y[n];
        } else {
          Vq = (1.0 - xtmp) * y[n] + xtmp * y[n + 1];
        }
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);

  return Vq;
}

//
// Arguments    : const double x[3]
// Return Type  : double
//
static double b_norm(const double x[3])
{
  double y;
  double scale;
  double absxk;
  double t;
  scale = 3.3121686421112381E-170;
  absxk = fabs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }

  absxk = fabs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = 1.0 + y * t * t;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  absxk = fabs(x[2]);
  if (absxk > scale) {
    t = scale / absxk;
    y = 1.0 + y * t * t;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  return scale * sqrtf(y);
}

//
// Arguments    : const double a[3]
//                double y[3]
// Return Type  : void
//
static void b_power(const double a[3], double y[3])
{
  y[0] = rt_powd_snf(a[0], 2.0);
  y[1] = rt_powd_snf(a[1], 2.0);
  y[2] = rt_powd_snf(a[2], 2.0);
}

//
// Arguments    : double x[4]
// Return Type  : void
//
static void b_sign(double x[4])
{
  double b_x;
  b_x = x[0];
  if (x[0] < 0.0) {
    b_x = -1.0;
  } else if (x[0] > 0.0) {
    b_x = 1.0;
  } else {
    if (x[0] == 0.0) {
      b_x = 0.0;
    }
  }

  x[0] = b_x;
  b_x = x[1];
  if (x[1] < 0.0) {
    b_x = -1.0;
  } else if (x[1] > 0.0) {
    b_x = 1.0;
  } else {
    if (x[1] == 0.0) {
      b_x = 0.0;
    }
  }

  x[1] = b_x;
  b_x = x[2];
  if (x[2] < 0.0) {
    b_x = -1.0;
  } else if (x[2] > 0.0) {
    b_x = 1.0;
  } else {
    if (x[2] == 0.0) {
      b_x = 0.0;
    }
  }

  x[2] = b_x;
  b_x = x[3];
  if (x[3] < 0.0) {
    b_x = -1.0;
  } else if (x[3] > 0.0) {
    b_x = 1.0;
  } else {
    if (x[3] == 0.0) {
      b_x = 0.0;
    }
  }

  x[3] = b_x;
}

//
// Arguments    : double *x
// Return Type  : void
//
static void b_sqrt(double *x)
{
  *x = sqrtf(*x);
}

//
// Arguments    : const double y[3]
//                const double x[3]
//                double r[3]
// Return Type  : void
//
static void c_atan2(const double y[3], const double x[3], double r[3])
{
  r[0] = rt_atan2d_snf(y[0], x[0]);
  r[1] = rt_atan2d_snf(y[1], x[1]);
  r[2] = rt_atan2d_snf(y[2], x[2]);
}

//
// Arguments    : const double varargin_1[2]
//                const double varargin_2[2]
//                double varargin_3
// Return Type  : double
//
static double c_interp1(const double varargin_1[2], const double varargin_2[2],
  double varargin_3)
{
  double Vq;
  double y_idx_0;
  double r;
  double y_idx_1;
  double x_idx_1;
  int k;
  int exitg1;
  y_idx_0 = varargin_2[0];
  r = varargin_1[0];
  y_idx_1 = varargin_2[1];
  x_idx_1 = varargin_1[1];
  k = 0;
  do {
    exitg1 = 0;
    if (k < 2) {
      if (isnan(varargin_1[k])) {
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

      Vq = NAN;
      if ((!isnan(varargin_3)) && (!(varargin_3 > x_idx_1)) && (!(varargin_3 <
            r))) {
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
// Arguments    : const double a[4]
//                double y[4]
// Return Type  : void
//
static void c_power(const double a[4], double y[4])
{
  y[0] = rt_powd_snf(a[0], 2.0);
  y[1] = rt_powd_snf(a[1], 2.0);
  y[2] = rt_powd_snf(a[2], 2.0);
  y[3] = rt_powd_snf(a[3], 2.0);
}

//
// Arguments    : double *x
// Return Type  : void
//
static void c_sign(double *x)
{
  if (*x < 0.0) {
    *x = -1.0;
  } else if (*x > 0.0) {
    *x = 1.0;
  } else {
    if (*x == 0.0) {
      *x = 0.0;
    }
  }
}

//
// Arguments    : double x[7]
// Return Type  : void
//
static void c_sqrt(double x[7])
{
  int k;
  for (k = 0; k < 7; k++) {
    x[k] = sqrtf(x[k]);
  }
}

//
// Arguments    : const double a[3]
//                const double b[3]
//                double c[3]
// Return Type  : void
//
static void cross(const double a[3], const double b[3], double c[3])
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

//
// Arguments    : const double y[4]
//                const double x[4]
//                double r[4]
// Return Type  : void
//
static void d_atan2(const double y[4], const double x[4], double r[4])
{
  r[0] = rt_atan2d_snf(y[0], x[0]);
  r[1] = rt_atan2d_snf(y[1], x[1]);
  r[2] = rt_atan2d_snf(y[2], x[2]);
  r[3] = rt_atan2d_snf(y[3], x[3]);
}

//
// Arguments    : double x[3]
// Return Type  : void
//
static void d_sqrt(double x[3])
{
  x[0] = sqrtf(x[0]);
  x[1] = sqrtf(x[1]);
  x[2] = sqrtf(x[2]);
}

//
// Arguments    : double x[4]
// Return Type  : void
//
static void e_sqrt(double x[4])
{
  x[0] = sqrtf(x[0]);
  x[1] = sqrtf(x[1]);
  x[2] = sqrtf(x[2]);
  x[3] = sqrtf(x[3]);
}

//
// Arguments    : const double varargin_1[10]
//                const double varargin_2[10]
//                double varargin_3
// Return Type  : double
//
static double interp1(const double varargin_1[10], const double varargin_2[10],
                      double varargin_3)
{
  double Vq;
  double y[10];
  double x[10];
  int low_i;
  int exitg1;
  double xtmp;
  int low_ip1;
  int high_i;
  int mid_i;
  memcpy(&y[0], &varargin_2[0], 10U * sizeof(double));
  memcpy(&x[0], &varargin_1[0], 10U * sizeof(double));
  low_i = 0;
  do {
    exitg1 = 0;
    if (low_i < 10) {
      if (isnan(varargin_1[low_i])) {
        exitg1 = 1;
      } else {
        low_i++;
      }
    } else {
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

      Vq = NAN;
      if ((!isnan(varargin_3)) && (!(varargin_3 > x[9])) && (!(varargin_3 < x
            [0]))) {
        low_i = 1;
        low_ip1 = 2;
        high_i = 10;
        while (high_i > low_ip1) {
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
        } else if (y[low_i - 1] == y[low_i]) {
          Vq = y[low_i - 1];
        } else {
          Vq = (1.0 - xtmp) * y[low_i - 1] + xtmp * y[low_i];
        }
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);

  return Vq;
}

//
// Arguments    : double V[10][11]
//                double Xq
//                double Yq
//                double extrapval
// Return Type  : double
//
static double interp2_dispatch(double V[10][11], double Xq, double Yq, double
  extrapval)
{
  double Vq;
  int low_i;
  int low_ip1;
  int high_i;
  int mid_i;
  static const double dv18[11] = { 0.0, 0.1, 0.2, 0.30000000000000004, 0.4, 0.5,
    0.6, 0.7, 0.8, 0.9, 1.0 };

  double dv19[10];
  double rx;
  double qx1;
  double d13;
  if ((Xq >= 0.0) && (Xq <= 1.0) && (Yq >= 0.0) && (Yq <= M_PI_2)) {
    low_i = 1;
    low_ip1 = 2;
    high_i = 11;
    while (high_i > low_ip1) {
      mid_i = (low_i + high_i) >> 1;
      if (Xq >= dv18[mid_i - 1]) {
        low_i = mid_i;
        low_ip1 = mid_i + 1;
      } else {
        high_i = mid_i;
      }
    }

    for (low_ip1 = 0; low_ip1 < 10; low_ip1++) {
      dv19[low_ip1] = 0.17453292519943295 * static_cast<double>(low_ip1);
    }

    low_ip1 = b_bsearch(dv19, Yq);
    rx = dv18[low_i - 1];
    if (Xq == rx) {
      qx1 = V[low_ip1 - 1][low_i - 1];
      Vq = V[low_ip1][low_i - 1];
    } else if (Xq == dv18[low_i]) {
      qx1 = V[low_ip1 - 1][low_i];
      Vq = V[low_ip1][low_i];
    } else {
      rx = (Xq - rx) / (dv18[low_i] - rx);
      if (V[low_ip1 - 1][low_i - 1] == V[low_ip1 - 1][low_i]) {
        qx1 = V[low_ip1 - 1][low_i - 1];
      } else {
        qx1 = (1.0 - rx) * V[low_ip1 - 1][low_i - 1] + rx * V[low_ip1 - 1][low_i];
      }

      if (V[low_ip1][low_i - 1] == V[low_ip1][low_i]) {
        Vq = V[low_ip1][low_i - 1];
      } else {
        Vq = (1.0 - rx) * V[low_ip1][low_i - 1] + rx * V[low_ip1][low_i];
      }
    }

    rx = 0.17453292519943295 * static_cast<double>((low_ip1 - 1));
    if ((Yq == rx) || (qx1 == Vq)) {
      Vq = qx1;
    } else {
      d13 = 0.17453292519943295 * static_cast<double>(low_ip1);
      if (!(Yq == d13)) {
        rx = (Yq - rx) / (d13 - rx);
        Vq = (1.0 - rx) * qx1 + rx * Vq;
      }
    }
  } else {
    Vq = extrapval;
  }

  return Vq;
}

//
// Arguments    : const double a[7]
//                double y[7]
// Return Type  : void
//
static void power(const double a[7], double y[7])
{
  int k;
  for (k = 0; k < 7; k++) {
    y[k] = rt_powd_snf(a[k], 2.0);
  }
}

//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  int b_u0;
  int b_u1;
  if (isnan(u0) || isnan(u1)) {
    y = NAN;
  } else if (isinf(u0) && isinf(u1)) {
    if (u0 > 0.0) {
      b_u0 = 1;
    } else {
      b_u0 = -1;
    }

    if (u1 > 0.0) {
      b_u1 = 1;
    } else {
      b_u1 = -1;
    }

    y = atan2((double)b_u0, (double)b_u1);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = M_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(M_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d15;
  double d16;
  if (isnan(u0) || isnan(u1)) {
    y = NAN;
  } else {
    d15 = fabs(u0);
    d16 = fabs(u1);
    if (isinf(u1)) {
      if (d15 == 1.0) {
        y = 1.0;
      } else if (d15 > 1.0) {
        if (u1 > 0.0) {
          y = INFINITY;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = INFINITY;
      }
    } else if (d16 == 0.0) {
      y = 1.0;
    } else if (d16 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrtf(u0);
    } else if ((u0 < 0.0) && (u1 > floorf(u1))) {
      y = NAN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

//
// Arguments    : double LAilDef
//                double ElevDef
//                double RudDef
//                double wIn
//                double u
//                double v
//                double w
//                double p
//                double q
//                double r
//                double *Fx
//                double *Fy
//                double *Fz
//                double *Mx
//                double *My
//                double *Mz
// Return Type  : void
//
void McFoamy_FM_v2(double LAilDef, double ElevDef, double RudDef, double wIn,
                   double u, double v, double w, double p, double q, double r,
                   double *Fx, double *Fy, double *Fz, double *Mx, double *My,
                   double *Mz)
{
  double B_rate[3];
  int i;
  int i0;
  double rods[32][7];
  static const short iv0[32][7] = { { -80, 48, 0, -255, 0, -80, 2 }, { -255, 0,
      -80, -328, 54, 0, 2 }, { -328, 54, 0, -425, 0, -85, 2 }, { -425, 0, -85,
      -540, 38, 0, 2 }, { -540, 38, 0, -645, 0, -60, 2 }, { -645, 0, -60, -695,
      115, 0, 2 }, { -745, 0, -93, -720, 115, 0, 2 }, { -80, 48, 0, -230, 0, 75,
      2 }, { -230, 0, 75, -185, 210, 0, 2 }, { -230, 0, 75, -265, 210, 0, 2 }, {
      -230, 0, 75, -330, 55, 0, 2 }, { -330, 55, 0, -445, 0, 65, 2 }, { -445, 0,
      65, -540, 35, 0, 2 }, { -540, 35, 0, -655, 0, 55, 2 }, { -655, 0, 55, -695,
      115, 0, 2 }, { -720, 115, 0, -745, 0, 65, 2 }, { -80, -48, 0, -255, 0, -80,
      2 }, { -255, 0, -80, -328, -54, 0, 2 }, { -328, -54, 0, -425, 0, -85, 2 },
    { -425, 0, -85, -540, -38, 0, 2 }, { -540, -38, 0, -645, 0, -60, 2 }, { -645,
      0, -60, -695, -115, 0, 2 }, { -745, 0, -93, -720, -115, 0, 2 }, { -80, -48,
      0, -230, 0, 75, 2 }, { -230, 0, 75, -185, -210, 0, 2 }, { -230, 0, 75,
      -265, -210, 0, 2 }, { -230, 0, 75, -330, -55, 0, 2 }, { -330, -55, 0, -445,
      0, 65, 2 }, { -445, 0, 65, -540, -35, 0, 2 }, { -540, -35, 0, -655, 0, 55,
      2 }, { -655, 0, 55, -695, -115, 0, 2 }, { -720, -115, 0, -745, 0, 65, 2 }
  };

  static const double dv0[32][6] = { { 160.0, 48.0, -6.61, -15.0, 0.0, -86.61 },
    { -15.0, 0.0, -86.61, -88.0, 54.0, -6.61 }, { -88.0, 54.0, -6.61, -185.0,
      0.0, -91.61 }, { -185.0, 0.0, -91.61, -300.0, 38.0, -6.61 }, { -300.0,
      38.0, -6.61, -405.0, 0.0, -66.61 }, { -405.0, 0.0, -66.61, -455.0, 115.0,
      -6.61 }, { -505.0, 0.0, -99.61, -480.0, 115.0, -6.61 }, { 160.0, 48.0,
      -6.61, 10.0, 0.0, 68.39 }, { 10.0, 0.0, 68.39, 55.0, 210.0, -6.61 }, {
      10.0, 0.0, 68.39, -25.0, 210.0, -6.61 }, { 10.0, 0.0, 68.39, -90.0, 55.0,
      -6.61 }, { -90.0, 55.0, -6.61, -205.0, 0.0, 58.39 }, { -205.0, 0.0, 58.39,
      -300.0, 35.0, -6.61 }, { -300.0, 35.0, -6.61, -415.0, 0.0, 48.39 }, { -
      415.0, 0.0, 48.39, -455.0, 115.0, -6.61 }, { -480.0, 115.0, -6.61, -505.0,
      0.0, 58.39 }, { 160.0, -48.0, -6.61, -15.0, 0.0, -86.61 }, { -15.0, 0.0,
      -86.61, -88.0, -54.0, -6.61 }, { -88.0, -54.0, -6.61, -185.0, 0.0, -91.61
    }, { -185.0, 0.0, -91.61, -300.0, -38.0, -6.61 }, { -300.0, -38.0, -6.61,
      -405.0, 0.0, -66.61 }, { -405.0, 0.0, -66.61, -455.0, -115.0, -6.61 }, { -
      505.0, 0.0, -99.61, -480.0, -115.0, -6.61 }, { 160.0, -48.0, -6.61, 10.0,
      0.0, 68.39 }, { 10.0, 0.0, 68.39, 55.0, -210.0, -6.61 }, { 10.0, 0.0,
      68.39, -25.0, -210.0, -6.61 }, { 10.0, 0.0, 68.39, -90.0, -55.0, -6.61 },
      { -90.0, -55.0, -6.61, -205.0, 0.0, 58.39 }, { -205.0, 0.0, 58.39, -300.0,
      -35.0, -6.61 }, { -300.0, -35.0, -6.61, -415.0, 0.0, 48.39 }, { -415.0,
      0.0, 48.39, -455.0, -115.0, -6.61 }, { -480.0, -115.0, -6.61, -505.0, 0.0,
      58.39 } };

  double vthr_x_tmp;
  double vthr_x;
  double vthr_y;
  double vthr_z;
  double d0;
  double d1;
  double d2;
  double psiT;
  double deltaT;
  double wOut;
  double CFx;
  double CFy;
  double CMy;
  double CMz;
  static const double dv1[10][11] = { { 0.156336454, 0.142939755, 0.12449467,
      0.103654093, 0.080437634, 0.054966543, 0.027444178, -0.001874328,
      -0.032728784, -0.064814458, -0.097502535 }, { 0.156336454, 0.143274298,
      0.125294853, 0.104978752, 0.082320505, 0.057426888, 0.03049664,
      0.001785855, -0.02842523, -0.059769732, -0.091725529 }, { 0.156336454,
      0.144236696, 0.127662714, 0.10890636, 0.087900516, 0.064729471,
      0.039588216, 0.012726044, -0.015539759, -0.044818343, -0.074619748 }, {
      0.156336454, 0.145715204, 0.131491476, 0.115283154, 0.096990697,
      0.076657686, 0.054473185, 0.030695378, 0.005645211, -0.020318541,
      -0.046776089 }, { 0.156336454, 0.147592844, 0.136533063, 0.123846679,
      0.109258684, 0.092816061, 0.074693054, 0.05514811, 0.034482474,
      0.013005062, -0.008957689 }, { 0.156336454, 0.149730239, 0.142227959,
      0.133924798, 0.124116893, 0.112598271, 0.099537509, 0.085260174,
      0.070037437, 0.054122603, 0.037754768 }, { 0.156336454, 0.151960378,
      0.148100411, 0.144287776, 0.139693524, 0.133992265, 0.127172755,
      0.119365045, 0.11072865, 0.101419857, 0.091641937 }, { 0.156336454,
      0.154120349, 0.153647038, 0.153875413, 0.153953862, 0.153533548,
      0.152570626, 0.151139537, 0.149362748, 0.147360112, 0.145232828 }, {
      0.156336454, 0.156099101, 0.158176177, 0.161071479, 0.163349211,
      0.165846339, 0.168638349, 0.171659606, 0.174874022, 0.17827724,
      0.181871525 }, { 0.156336454, 0.157836404, 0.161723393, 0.163852071,
      0.165797264, 0.168474477, 0.171699702, 0.175259993, 0.179027738,
      0.182929631, 0.18696977 } };

  static const double dv2[10][11] = { { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0 }, { 0.0, 0.000304035, 0.000672049, 0.001105525, 0.001600104,
      0.00214306, 0.002695169, 0.003243617, 0.003737499, 0.004001698,
      0.003755281 }, { 0.0, 0.000593774, 0.001308186, 0.002139399, 0.003086457,
      0.004124254, 0.005206331, 0.0062742, 0.007197636, 0.007738064, 0.007554928
    }, { 0.0, 0.000855707, 0.001875984, 0.003046099, 0.004372415, 0.005829393,
      0.007370573, 0.008920992, 0.010310416, 0.011301405, 0.011576152 }, { 0.0,
      0.001081293, 0.002343891, 0.003777072, 0.005383801, 0.007152135,
      0.009046337, 0.010988567, 0.012837704, 0.014405295, 0.015455208 }, { 0.0,
      0.001265997, 0.002685232, 0.004271303, 0.006042688, 0.007991019,
      0.010092937, 0.012301024, 0.014518175, 0.016623841, 0.018472368 }, { 0.0,
      0.001410344, 0.002915524, 0.004529052, 0.006286839, 0.008217389,
      0.010327588, 0.012588497, 0.014938663, 0.017317246, 0.019651635 }, { 0.0,
      0.001515466, 0.003073033, 0.004674445, 0.006333996, 0.008079318,
      0.009939438, 0.01193205, 0.014050094, 0.016267158, 0.018549771 }, { 0.0,
      0.001581189, 0.003203695, 0.004870829, 0.006596278, 0.008340642,
      0.010095901, 0.011869291, 0.013673189, 0.015518142, 0.017413651 }, { 0.0,
      0.00160638, 0.003296003, 0.005159087, 0.007149694, 0.00921631, 0.01134185,
      0.013520997, 0.015752967, 0.018038805, 0.020375896 } };

  double dv3[10];
  static const double dv4[10] = { -0.097502535, -0.091725529, -0.074619748,
    -0.046776089, -0.008957689, 0.037754768, 0.091641937, 0.145232828,
    0.181871525, 0.18696977 };

  static const double dv5[10][11] = { { -0.009379255, -0.009414602, -0.00922517,
      -0.008715215, -0.007778011, -0.006305113, -0.004198032, -0.001374064,
      0.002229593, 0.006644511, 0.011790387 }, { -0.009379255, -0.009414258,
      -0.009234779, -0.008749876, -0.007856254, -0.0064465, -0.004423154,
      -0.001706255, 0.001762192, 0.00599214, 0.01091461 }, { -0.009379255,
      -0.009412543, -0.009259831, -0.008846651, -0.008075885, -0.006846983,
      -0.005068319, -0.002664393, 0.00040563, 0.004131457, 0.008431755 }, { -
      0.009379255, -0.009408081, -0.009292432, -0.008982547, -0.008395035,
      -0.007438956, -0.006031948, -0.004110867, -0.001648602, 0.001338374,
      0.004772847 }, { -0.009379255, -0.009399845, -0.0093192, -0.009121875,
      -0.008745071, -0.008111231, -0.007148975, -0.005807856, -0.004067627,
      -0.001941328, 0.000516475 }, { -0.009379255, -0.009386461, -0.009327092,
      -0.009224484, -0.009041916, -0.008720706, -0.008203888, -0.007449251,
      -0.006438264, -0.005174359, -0.003684646 }, { -0.009379255, -0.009367587,
      -0.009306845, -0.009256742, -0.009209727, -0.009127023, -0.008967104,
      -0.00869943, -0.008307869, -0.007785904, -0.007140396 }, { -0.009379255,
      -0.009343565, -0.009256742, -0.009206982, -0.009213159, -0.00925331,
      -0.009300668, -0.009331211, -0.009330867, -0.009294491, -0.009222768 }, {
      -0.009379255, -0.009317827, -0.009186392, -0.009150015, -0.00927596,
      -0.009427299, -0.009571088, -0.009705269, -0.009834645, -0.009963678,
      -0.010097172 }, { -0.009379255, -0.009292775, -0.009135945, -0.009359351,
      -0.009773903, -0.01020081, -0.010605068, -0.010994226, -0.011382011,
      -0.011781807, -0.012199449 } };

  double a;
  static const double dv6[10][11] = { { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0 }, { 0.0, 0.000624918, 0.00118875, 0.001686008, 0.00210811,
      0.002447508, 0.002695622, 0.002844902, 0.002886426, 0.00279274,
      0.002507907 }, { 0.0, 0.001228902, 0.002355195, 0.003352112, 0.004209014,
      0.004910802, 0.005441004, 0.005784863, 0.005919044, 0.005797904,
      0.005340455 }, { 0.0, 0.001786557, 0.003475311, 0.004976005, 0.006291729,
      0.007398462, 0.008275611, 0.008902245, 0.009246104, 0.009258801,
      0.008857975 }, { 0.0, 0.002280383, 0.004508604, 0.006530234, 0.008331891,
      0.009900534, 0.011210425, 0.012237884, 0.012952713, 0.013314074,
      0.013262255 }, { 0.0, 0.002695622, 0.005373056, 0.007916309, 0.010270131,
      0.012377212, 0.014212843, 0.015771877, 0.017029948, 0.017960633,
      0.018530643 }, { 0.0, 0.003025754, 0.006031948, 0.008970193, 0.01181784,
      0.01452513, 0.017039214, 0.019318224, 0.021326814, 0.023040962,
      0.024451059 }, { 0.0, 0.003265632, 0.006465718, 0.009615701, 0.012740975,
      0.015837766, 0.018875875, 0.021818238, 0.024630882, 0.027285666,
      0.029763374 }, { 0.0, 0.003407706, 0.006690153, 0.009827095, 0.012792794,
      0.015769475, 0.018797288, 0.021878979, 0.025002195, 0.028151834,
      0.031305935 }, { 0.0, 0.003447857, 0.006735452, 0.009570745, 0.012191899,
      0.014745447, 0.017259188, 0.019738611, 0.022188521, 0.024615439,
      0.027038239 } };

  static const double dv7[10][11] = { { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0 }, { 0.0, 0.000466715, 0.000650999, 0.000651342, 0.000542556,
      0.000378176, 0.000193549, 7.5498E-6, -0.00017193, -0.00033528, -0.00047049
    }, { 0.0, 0.000928969, 0.001316754, 0.001337687, 0.001141049, 0.000833567,
      0.00048456, 0.000133151, -0.000201786, -0.000507209, -0.000764589 }, { 0.0,
      0.001379898, 0.002009277, 0.002091295, 0.001847642, 0.001435835,
      0.000959511, 0.000477353, 1.9904E-5, -0.000396021, -0.000754637 }, { 0.0,
      0.001812295, 0.002733028, 0.002939961, 0.002706947, 0.002244007,
      0.001687037, 0.001115998, 0.000571383, 7.41253E-5, -0.000362734 }, { 0.0,
      0.002217926, 0.003466731, 0.003887461, 0.003752251, 0.003306126,
      0.002723076, 0.002106738, 0.001510303, 0.000958825, 0.000464313 }, { 0.0,
      0.002588209, 0.004186021, 0.004884035, 0.004944433, 0.004618762,
      0.004095081, 0.003494871, 0.002888485, 0.002311955, 0.001782782 }, { 0.0,
      0.002912507, 0.004849031, 0.005825014, 0.006105387, 0.005966058,
      0.005597491, 0.005116363, 0.004593367, 0.004069343, 0.003567967 }, { 0.0,
      0.003173662, 0.005419041, 0.00627663, 0.005895365, 0.005365163,
      0.004922813, 0.004597142, 0.004365844, 0.004200435, 0.004069686 }, { 0.0,
      0.003355543, 0.005744369, 0.005094056, 0.003495558, 0.001968096,
      0.000671589, -0.000392246, -0.001249149, -0.001923827, -0.002446135 } };

  static const double dv8[10] = { 0.0, 0.003755281, 0.007554928, 0.011576152,
    0.015455208, 0.018472368, 0.019651635, 0.018549771, 0.017413651, 0.020375896
  };

  double b_a;
  double MX;
  static const double dv9[10] = { 0.011790387, 0.01091461, 0.008431755,
    0.004772847, 0.000516475, -0.003684646, -0.007140396, -0.009222768,
    -0.010097172, -0.012199449 };

  double c_a;
  static const double dv10[10] = { 0.0, 0.002507907, 0.005340455, 0.008857975,
    0.013262255, 0.018530643, 0.024451059, 0.029763374, 0.031305935, 0.027038239
  };

  static const double dv11[10] = { 0.0, -0.00047049, -0.000764589, -0.000754637,
    -0.000362734, 0.000464313, 0.001782782, 0.003567967, 0.004069686,
    -0.002446135 };

  double ThrForce_idx_0;
  double ThrForce_idx_1_tmp;
  double b_ThrForce_idx_1_tmp;
  double ThrForce_idx_1;
  double ThrForce_idx_2;
  static const double dv12[50] = { 0.21706999999999999, 0.21139, 0.21543,
    0.21218, 0.20885, 0.20550000000000002, 0.20226, 0.71127, 0.71972, 0.71515,
    0.73394, 0.73581, 0.74414, 0.76076, 0.20509, 0.20509, 0.20509, 0.20509,
    0.1675, 0.2915, 0.3765, 0.4825, 0.5925, 0.67, 0.7325, 0.155,
    0.20750000000000002, 0.2475, 0.28, 0.3875, 0.4925, 0.5975, 0.675, 0.7325,
    0.1675, 0.2915, 0.3765, 0.4825, 0.5925, 0.67, 0.7325, 0.155,
    0.20750000000000002, 0.2475, 0.28, 0.3875, 0.4925, 0.5975, 0.675, 0.7325 };

  static const double dv13[50] = { 0.03338, 0.091760000000000008, 0.15114,
    0.21201, 0.27455, 0.33742, 0.40021, 0.03953, 0.09333, 0.15156999999999998,
    0.02354, 0.02996, 0.094629999999999992, 0.14856, 0.053810000000000004,
    0.01794, 0.01863, 0.05589, 0.046647615158762409, 0.048259714048054618,
    0.050351266121121521, 0.0465537323960174, 0.035510561809129405,
    0.064855608855364233, 0.0739493069609175, 0.044522466238967488,
    0.1114955156048888, 0.1114955156048888, 0.046502688094345684,
    0.042573465914816006, 0.0369120576505835, 0.032596012026013241,
    0.063737743919909814, 0.066049224068114532, 0.046647615158762409,
    0.048259714048054618, 0.050351266121121521, 0.0465537323960174,
    0.035510561809129405, 0.064855608855364233, 0.0739493069609175,
    0.044522466238967488, 0.1114955156048888, 0.1114955156048888,
    0.046502688094345684, 0.042573465914816006, 0.0369120576505835,
    0.032596012026013241, 0.063737743919909814, 0.066049224068114532 };

  double VProp_Axial[50];
  double vts_xcomp[3];
  double vws_xcomp[7];
  static const double b[7] = { 0.03338, 0.091760000000000008, 0.15114, 0.21201,
    0.27455, 0.33742, 0.40021 };

  static const double b_b[7] = { 0.022930000000000006, 0.028609999999999997,
    0.024569999999999981, 0.027819999999999984, 0.031149999999999983,
    0.034499999999999975, 0.037739999999999996 };

  double vts_zcomp[3];
  double vws_zcomp[7];
  double vtp_xcomp[3];
  double vwp_xcomp[7];
  static const double c_b[7] = { -0.03338, -0.091760000000000008, -0.15114,
    -0.21201, -0.27455, -0.33742, -0.40021 };

  double vtp_zcomp[3];
  double vwp_zcomp[7];
  double Mx_r[4];
  double vr_xcomp[4];
  double vr_ycomp[4];
  double vB_xcomp[4];
  double vB_ycomp[4];
  double a_ws[7];
  double a_wp[7];
  double a_ts[3];
  double a_tp[3];
  double a_r[4];
  double a_B[4];
  static const double dv14[7] = { 0.0, 0.10204, 0.099230000000000013, 0.09636,
    0.093409999999999993, 0.09044, 0.08749 };

  static const double dv15[7] = { 0.2534, 0.23758, 0.22269, 0.20742,
    0.19174000000000002, 0.17597000000000002, 0.16028 };

  double d3;
  double CL_ws[7];
  double CD_ws[7];
  double CM_ws[7];
  double CL_ts[3];
  double CL_wp[7];
  double CD_wp[7];
  double d4;
  double CM_wp[7];
  double CL_tp[3];
  double CD_tp_idx_0;
  double CM_tp_idx_0;
  double CD_ts_idx_1;
  double CM_ts_idx_1;
  double CD_tp_idx_1;
  double CM_tp_idx_1;
  double d5;
  double d6;
  double d7;
  double CL_r[4];
  double CD_r_idx_0;
  double CM_r_idx_0;
  double d8;
  double d9;
  double d10;
  double CL_B[4];
  double CD_B_idx_0;
  double CM_B_idx_0;
  double CD_r_idx_1;
  double CM_r_idx_1;
  double CD_B_idx_1;
  double CM_B_idx_1;
  double CD_r_idx_2;
  double CM_r_idx_2;
  double CD_B_idx_2;
  double CM_B_idx_2;
  double D_rods_idx_0;
  double M_rods_idx_0;
  double D_rods_idx_1;
  double M_rods_idx_1;
  double D_rods_idx_2;
  double M_rods_idx_2;
  double dv16[7];
  double ri[3];
  double Mz_ws[7];
  double li[3];
  double My_ws[7];
  double Vi[3];
  static const double Fx_ws_tmp[7] = { 0.00888718145, 0.0085448022800000013,
    0.008087015186250001, 0.007950460455, 0.0073494421350000011,
    0.006817187781250001, 0.0060768158500000011 };

  double V_l_l[3];
  double Fx_ws[7];
  double Mz_ts[3];
  double Mx_wp[7];
  double dv17[4];
  double CtrlDef[4];
  double d11;
  double d12;
  static const double Cw[7] = { 0.2534, 0.23758, 0.22269, 0.20742,
    0.19174000000000002, 0.17597000000000002, 0.16028 };

  B_rate[0] = p;
  B_rate[1] = q;
  B_rate[2] = r;

  //  MCFOAMY GEOMETRY
  //         x     y    z
  //  CG = [-238.58, 0, 5.89]; %measured from propeller plane using CAD
  //  CG = [-270, 0, 5.89]; % CG_x measured roughly by hand (picking plane up)
  //  CG_x measured by hand, CG_z from CAD
  //          #segs  wing root   wing tip  wing span
  //  Wing Sections
  //             Span    Area     Chord FlapChord   x     y     z
  //  Tail Sections
  //             Span    Area     Chord FlapChord   x     y     z
  //  Rudder Sections
  //             Span    Area     Chord FlapChord   x     y     z
  //  Body Sections
  //             Span    Area    Chord  FlapChord   x     y     z
  //  reinforcement rods
  //  small rod diameter
  //  structural rods (7x2mm top)
  //  9x2mm rods bottom
  //  mirror it for y
  //  add the diameter in last row
  //  rod centers [mm] w/r to prop plane
  //  shift it by CG
  for (i = 0; i < 32; i++) {
    for (i0 = 0; i0 < 7; i0++) {
      rods[i][i0] = iv0[i][i0];
    }

    for (i0 = 0; i0 < 6; i0++) {
      rods[i][i0] = dv0[i][i0];
    }

    for (i0 = 0; i0 < 7; i0++) {
      rods[i][i0] *= 0.001;
    }
  }

  //  convert to [m]
  //  Geometric properties of components
  //  AIRCRAFT CG
  //  Measured w.r.t the nose (from propeller plane)
  //  THRUSTER
  //  Thrust is assumed to act in x direction only
  //  WING GEOMETRY
  //  Aspect ratio
  //  TAIL GEOMETRY
  //  Aspect ratio
  //  RUDDER GEOMETRY
  //  Aspect ratio
  //  BODY GEOMETRY
  //  Aspect ratio
  //  THRUSTER MODEL
  //  PROPELLER MAPS FOR YAK54 (ELECTRIFLY 10X4.7 PROP)
  //  from thruster model
  //  Negative J data from Selig's database (for APC 10x4.7) % we don't use it
  //  CALCULATIONS
  //  negative b/c it is from c.g. to propeller plane
  vthr_x_tmp = u + q * -0.00661;
  vthr_x = vthr_x_tmp - r * -0.0;
  vthr_y = (v + r * 0.24) - p * -0.00661;
  vthr_z = (w + p * -0.0) - q * 0.24;

  // -------------------------------------------------------------
  //  total velocity
  //  in-plane velocity
  d0 = vthr_y * vthr_y;
  d1 = vthr_z * vthr_z;
  d2 = d0 + d1;
  b_sqrt(&d2);
  psiT = rt_atan2d_snf(d2, fabs(vthr_x));

  //  azimuth angle of thruster (0 to +90 deg.)
  deltaT = rt_atan2d_snf(vthr_z, vthr_y);

  //  in-plane angle (0 to +/- 180)
  if (wIn < 1716.0) {
    //  Too slow to run the motor
    wOut = 0.0;
    CFx = 0.0;
    CFy = 0.0;
    vthr_z = 0.0;
    CMy = 0.0;
    CMz = 0.0;
  } else {
    wOut = wIn;
    d0 = (vthr_x * vthr_x + d0) + d1;
    b_sqrt(&d0);
    vthr_y = d0 / (wIn / 60.0 * 0.254);

    //  advance ratio based on total velocity
    if (vthr_x >= 0.0) {
      //  Forward flight
      if (vthr_y <= 1.0) {
        CFx = interp2_dispatch(*(double (*)[10][11])&dv1[0][0], vthr_y, psiT,
          NAN);
        CFy = interp2_dispatch(*(double (*)[10][11])&dv2[0][0], vthr_y, psiT,
          NAN);
        vthr_z = interp2_dispatch(*(double (*)[10][11])&dv5[0][0], vthr_y, psiT,
          NAN);
        CMy = interp2_dispatch(*(double (*)[10][11])&dv6[0][0], vthr_y, psiT,
          NAN);
        CMz = interp2_dispatch(*(double (*)[10][11])&dv7[0][0], vthr_y, psiT,
          NAN);
      } else {
        for (i = 0; i < 10; i++) {
          dv3[i] = 0.17453292519943295 * static_cast<double>(i);
        }

        CFx = interp1(dv3, dv4, psiT);
        for (i = 0; i < 10; i++) {
          dv3[i] = 0.17453292519943295 * static_cast<double>(i);
        }

        CFy = interp1(dv3, dv8, psiT);
        for (i = 0; i < 10; i++) {
          dv3[i] = 0.17453292519943295 * static_cast<double>(i);
        }

        vthr_z = interp1(dv3, dv9, psiT);
        for (i = 0; i < 10; i++) {
          dv3[i] = 0.17453292519943295 * static_cast<double>(i);
        }

        CMy = interp1(dv3, dv10, psiT);
        for (i = 0; i < 10; i++) {
          dv3[i] = 0.17453292519943295 * static_cast<double>(i);
        }

        CMz = interp1(dv3, dv11, psiT);
      }
    } else {
      //  Rearward flight
      //  A/c to Bart's paper Fig.8 & 9, for rearward flight, Cfx and Cmx have
      //  approx. static value
      CFx = 0.15633;
      CFy = 0.0;
      vthr_z = -0.00938;
      CMy = interp2_dispatch(*(double (*)[10][11])&dv6[0][0], vthr_y, psiT,
        NAN);
      CMz = interp2_dispatch(*(double (*)[10][11])&dv7[0][0], vthr_y, psiT,
        NAN);
    }
  }

  //  Transformation into UAV XYZ frame
  //  Aerodynamic Forces and Moments
  a = wOut / 60.0;
  psiT = wOut / 60.0;
  vthr_y = wOut / 60.0;
  b_a = wOut / 60.0;
  MX = 1.225 * (b_a * b_a) * 0.001057227821024 * vthr_z;
  b_a = wOut / 60.0;
  c_a = wOut / 60.0;

  //  Induced Velocity at Propeller Plane
  //  Induced velocity in hover
  d0 = 0.15633;
  b_sqrt(&d0);
  d1 = 1.59 * (wOut / 60.0) * 0.254;
  if ((vthr_x >= -0.2 * (0.5 * (d1 * d0))) && (CFx > 0.0)) {
    //  if reverse velocity is greater than 20% of hover velocity
    d0 = CFx;
    b_sqrt(&d0);
    vthr_z = 0.5 * (d1 * d0);
    MX *= 0.2;
  } else {
    vthr_z = 0.0;
  }

  //  Swirling Flow
  MX *= 0.4;
  ThrForce_idx_0 = 1.225 * (a * a) * 0.004162314256 * CFx;
  ThrForce_idx_1_tmp = sinf(deltaT);
  b_ThrForce_idx_1_tmp = cosf(deltaT);
  ThrForce_idx_1 = 0.5 * (1.225 * (psiT * psiT) * 0.004162314256 * (-CFy *
    b_ThrForce_idx_1_tmp + 0.0 * ThrForce_idx_1_tmp));
  ThrForce_idx_2 = 0.5 * (1.225 * (vthr_y * vthr_y) * 0.004162314256 * (-CFy *
    ThrForce_idx_1_tmp - 0.0 * b_ThrForce_idx_1_tmp));

  //  GYROSCOPIC MOMENTS
  //  Rotor moment of inertia
  //  PROPELLER SLIPSTREAM MODEL
  //  Equations used from propeller slipstream paper
  //  Eq. 11
  //  Eq. 12
  vthr_y = 2.0 * vthr_z * 0.91823899371069173;

  //  Efflux velocity from Vi0_avg
  for (i = 0; i < 50; i++) {
    if (dv12[i] < 0.194056) {
      a = dv12[i] / 0.127;
      d0 = 1.0 + a * a;
      b_sqrt(&d0);
      d0 = vthr_z * (1.0 + dv12[i] / 0.127 / d0);
      VProp_Axial[i] = d0;

      //  Eq. 1
    } else if (dv12[i] < 0.319532) {
      //  Eq. 13
      //  Eq. 14
      a = (dv13[i] - 0.058946599999999995 * (1.0 - 0.1294 * (dv12[i] - 0.194056)
            / 0.18796)) / (0.05210289974 + 0.1326 * ((dv12[i] - 0.194056) -
        0.09398));
      d0 = vthr_y * (1.24 - 0.0765 * (dv12[i] - 0.194056) / 0.18796) * expf(
        -(a * a));
      VProp_Axial[i] = d0;

      //  Eq. 15
    } else {
      //  Eq. 16
      //  Eq. 17
      a = (dv13[i] - 0.058946599999999995 * (1.3 - 0.3059 * (dv12[i] - 0.194056)
            / 0.18796)) / (0.030510760159999994 + 0.2295 * ((dv12[i] - 0.194056)
        - 0.09398));
      d0 = vthr_y * (1.37 - 0.1529 * (dv12[i] - 0.194056) / 0.18796) * expf(
        -(a * a));
      VProp_Axial[i] = d0;

      //  Eq. 18
    }

    if (d0 < 0.01) {
      VProp_Axial[i] = 0.0;
    }
  }

  //  AERODYNAMICS
  // --------------------------------------------------------------------------
  //                        Aerodynamic Constants
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------
  //                        Wind disturbances
  // --------------------------------------------------------------------------
  //  Wind disturbance in body frame
  // --------------------------------------------------------------------------
  //      Body axis components of vehicle velocity relative to the air
  // --------------------------------------------------------------------------
  //  Assuming wind velocity is in body frame
  // --------------------------------------------------------------------------
  //               Slipstream velocities at reference points
  // --------------------------------------------------------------------------
  //  rods
  // --------------------------------------------------------------------------
  //                              Control Inputs
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------
  //                    Wing, Tail, Rudder and Body velocities
  // --------------------------------------------------------------------------
  for (i = 0; i < 7; i++) {
    vws_xcomp[i] = (vthr_x_tmp - r * b[i]) + VProp_Axial[i];
    d0 = q * b_b[i];
    vws_zcomp[i] = (w + p * b[i]) - d0;
    vwp_xcomp[i] = (vthr_x_tmp - r * c_b[i]) + VProp_Axial[i];
    d0 = (w + p * c_b[i]) - d0;
    vwp_zcomp[i] = d0;
  }

  vts_xcomp[0] = (vthr_x_tmp - r * 0.03953) + VProp_Axial[7];
  d0 = q * -0.47126999999999997;
  vts_zcomp[0] = (w + p * 0.03953) - d0;
  vtp_xcomp[0] = (vthr_x_tmp - r * -0.03953) + VProp_Axial[7];
  d0 = (w + p * -0.03953) - d0;
  vtp_zcomp[0] = d0;
  vts_xcomp[1] = (vthr_x_tmp - r * 0.09333) + VProp_Axial[8];
  d0 = q * -0.47972000000000004;
  vts_zcomp[1] = (w + p * 0.09333) - d0;
  vtp_xcomp[1] = (vthr_x_tmp - r * -0.09333) + VProp_Axial[8];
  d0 = (w + p * -0.09333) - d0;
  vtp_zcomp[1] = d0;
  vts_xcomp[2] = (vthr_x_tmp - r * 0.15156999999999998) + VProp_Axial[9];
  d0 = q * -0.47514999999999996;
  vts_zcomp[2] = (w + p * 0.15156999999999998) - d0;
  vtp_xcomp[2] = (vthr_x_tmp - r * -0.15156999999999998) + VProp_Axial[9];
  d0 = (w + p * -0.15156999999999998) - d0;
  vtp_zcomp[2] = d0;
  d0 = r * 0.0;
  Mx_r[0] = d0;
  vr_xcomp[0] = ((u + q * 0.016929999999999997) - d0) + VProp_Axial[10];
  vr_ycomp[0] = 0.02354;
  d0 = r * 0.0;
  Mx_r[1] = d0;
  vr_xcomp[1] = ((u + q * -0.03657) - d0) + VProp_Axial[11];
  vr_ycomp[1] = -0.029959999999999997;
  d0 = r * 0.0;
  Mx_r[2] = d0;
  vr_xcomp[2] = ((u + q * -0.10124) - d0) + VProp_Axial[12];
  vr_ycomp[2] = -0.094629999999999992;
  d0 = r * 0.0;
  vr_xcomp[3] = ((u + q * -0.15517) - d0) + VProp_Axial[13];
  vr_ycomp[3] = -0.14856;
  b_sign(vr_ycomp);
  vr_ycomp[0] = ((v + r * -0.49394000000000005) - p * 0.016929999999999997) +
    vr_ycomp[0] * 0.0;
  vB_xcomp[0] = ((u + q * 0.047200000000000006) - Mx_r[0]) + VProp_Axial[14];
  vB_ycomp[0] = 0.053810000000000004;
  vr_ycomp[1] = ((v + r * -0.49581) - p * -0.03657) + vr_ycomp[1] * 0.0;
  vB_xcomp[1] = ((u + q * 0.01133) - Mx_r[1]) + VProp_Axial[15];
  vB_ycomp[1] = 0.01794;
  vr_ycomp[2] = ((v + r * -0.50414) - p * -0.10124) + vr_ycomp[2] * 0.0;
  vB_xcomp[2] = ((u + q * -0.025240000000000002) - Mx_r[2]) + VProp_Axial[16];
  vB_ycomp[2] = -0.01863;
  vr_ycomp[3] = ((v + r * -0.52076) - p * -0.15517) + vr_ycomp[3] * 0.0;
  vB_xcomp[3] = ((u + q * -0.0625) - d0) + VProp_Axial[17];
  vB_ycomp[3] = -0.05589;
  b_sign(vB_ycomp);
  vthr_y = v + r * 0.03491;
  vB_ycomp[0] = (vthr_y - p * 0.047200000000000006) + vB_ycomp[0] * 0.0;
  vB_ycomp[1] = (vthr_y - p * 0.01133) + vB_ycomp[1] * 0.0;
  vB_ycomp[2] = (vthr_y - p * -0.025240000000000002) + vB_ycomp[2] * 0.0;
  vB_ycomp[3] = (vthr_y - p * -0.0625) + vB_ycomp[3] * 0.0;

  // --------------------------------------------------------------------------
  //               Wing, Tail, Rudder and Body angles of attack
  // --------------------------------------------------------------------------
  //  Angle of attack = atan(V_zcomp / V_xcomp)
  //  range is -180 -> 180
  b_atan2(vws_zcomp, vws_xcomp, a_ws);
  b_atan2(vwp_zcomp, vwp_xcomp, a_wp);
  c_atan2(vts_zcomp, vts_xcomp, a_ts);
  c_atan2(vtp_zcomp, vtp_xcomp, a_tp);

  //  Vertical AoA for rudder and body = atan(V_ycomp / V_xcomp)
  //  range is -180 -> 180
  d_atan2(vr_ycomp, vr_xcomp, a_r);
  d_atan2(vB_ycomp, vB_xcomp, a_B);

  // --------------------------------------------------------------------------
  //                    Wing lift and drag coefficient
  // --------------------------------------------------------------------------
  for (i = 0; i < 7; i++) {
    McFoamy_Airfoil_Simplified(a_ws[i], dv14[i], dv15[i], -LAilDef,
      2.0920035851844676, 0.02, 1.98, 0.0, &d0, &d1, &d2, &d3);
    CL_ws[i] = d1;
    CD_ws[i] = d2;
    CM_ws[i] = d3;
    McFoamy_Airfoil_Simplified(a_wp[i], dv14[i], dv15[i], LAilDef,
      2.0920035851844676, 0.02, 1.98, 0.0, &d0, &d1, &d2, &d3);
    CL_wp[i] = d1;
    CD_wp[i] = d2;
    CM_wp[i] = d3;
  }

  // --------------------------------------------------------------------------
  //                   Tail's lift and drag coefficient
  // --------------------------------------------------------------------------
  McFoamy_Airfoil_Simplified(a_ts[0], 0.081870000000000012, 0.11741, ElevDef,
    1.2420750265020692, 0.02, 1.98, 0.0, &d0, &d1, &d2, &d3);
  CL_ts[0] = d1;
  vthr_x = d2;
  deltaT = d3;
  McFoamy_Airfoil_Simplified(a_tp[0], 0.081870000000000012, 0.11741, ElevDef,
    1.2420750265020692, 0.02, 1.98, 0.0, &d0, &CFy, &vthr_x_tmp, &d4);
  CL_tp[0] = CFy;
  CD_tp_idx_0 = vthr_x_tmp;
  CM_tp_idx_0 = d4;
  McFoamy_Airfoil_Simplified(a_ts[1], 0.11148999999999999, 0.14564, ElevDef,
    1.2420750265020692, 0.02, 1.98, 0.0, &d0, &d1, &d2, &d3);
  CL_ts[1] = d1;
  CD_ts_idx_1 = d2;
  CM_ts_idx_1 = d3;
  McFoamy_Airfoil_Simplified(a_tp[1], 0.11148999999999999, 0.14564, ElevDef,
    1.2420750265020692, 0.02, 1.98, 0.0, &d0, &CFy, &vthr_x_tmp, &d4);
  CL_tp[1] = CFy;
  CD_tp_idx_1 = vthr_x_tmp;
  CM_tp_idx_1 = d4;
  McFoamy_Airfoil_Simplified(a_ts[2], 0.14316, 0.14316, ElevDef,
    1.2420750265020692, 0.02, 1.98, 0.0, &d0, &d1, &d2, &d3);
  McFoamy_Airfoil_Simplified(a_tp[2], 0.14316, 0.14316, ElevDef,
    1.2420750265020692, 0.02, 1.98, 0.0, &d0, &CFy, &vthr_x_tmp, &d4);

  // --------------------------------------------------------------------------
  //                   Rudder lift and drag coefficient
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------
  //                   Fuselage lift and drag coefficient
  // --------------------------------------------------------------------------
  //  check sign of deflection into flapped code
  McFoamy_Airfoil_Simplified(a_r[0], 0.12036000000000001, 0.20062000000000002,
    -RudDef, 1.2984123165744528, 0.02, 1.98, 0.0, &d0, &d5, &d6, &d7);
  CL_r[0] = d5;
  CD_r_idx_0 = d6;
  CM_r_idx_0 = d7;
  McFoamy_Airfoil_Simplified(a_B[0], 0.0, 0.64171, 0.0, 0.22793785354755264,
    0.02, 1.98, 0.0, &d0, &d8, &d9, &d10);
  CL_B[0] = d8;
  CD_B_idx_0 = d9;
  CM_B_idx_0 = d10;

  //  check sign of deflection into flapped code
  McFoamy_Airfoil_Simplified(a_r[1], 0.11269, 0.19294999999999998, -RudDef,
    1.2984123165744528, 0.02, 1.98, 0.0, &d0, &d5, &d6, &d7);
  CL_r[1] = d5;
  CD_r_idx_1 = d6;
  CM_r_idx_1 = d7;
  McFoamy_Airfoil_Simplified(a_B[1], 0.0, 0.64171, 0.0, 0.22793785354755264,
    0.02, 1.98, 0.0, &d0, &d8, &d9, &d10);
  CL_B[1] = d8;
  CD_B_idx_1 = d9;
  CM_B_idx_1 = d10;

  //  check sign of deflection into flapped code
  McFoamy_Airfoil_Simplified(a_r[2], 0.10343000000000001, 0.16949, -RudDef,
    1.2984123165744528, 0.02, 1.98, 0.0, &d0, &d5, &d6, &d7);
  CL_r[2] = d5;
  CD_r_idx_2 = d6;
  CM_r_idx_2 = d7;
  McFoamy_Airfoil_Simplified(a_B[2], 0.0, 0.64171, 0.0, 0.22793785354755264,
    0.02, 1.98, 0.0, &d0, &d8, &d9, &d10);
  CL_B[2] = d8;
  CD_B_idx_2 = d9;
  CM_B_idx_2 = d10;

  //  check sign of deflection into flapped code
  McFoamy_Airfoil_Simplified(a_r[3], 0.14114, 0.14114, -RudDef,
    1.2984123165744528, 0.02, 1.98, 0.0, &d0, &d5, &d6, &d7);
  McFoamy_Airfoil_Simplified(a_B[3], 0.0, 0.64171, 0.0, 0.22793785354755264,
    0.02, 1.98, 0.0, &d0, &d8, &d9, &d10);

  // --------------------------------------------------------------------------
  //                            Structural rods
  // --------------------------------------------------------------------------
  D_rods_idx_0 = 0.0;
  M_rods_idx_0 = 0.0;
  D_rods_idx_1 = 0.0;
  M_rods_idx_1 = 0.0;
  D_rods_idx_2 = 0.0;
  M_rods_idx_2 = 0.0;
  for (i = 0; i < 32; i++) {
    //  pt1
    //  pt2
    //  diameter of the rod
    //  center of the rod w/r to the CG
    ri[0] = (rods[i][0] + rods[i][3]) / 2.0;
    li[0] = rods[i][3] - rods[i][0];
    ri[1] = (rods[i][1] + rods[i][4]) / 2.0;
    li[1] = rods[i][4] - rods[i][1];
    ri[2] = (rods[i][2] + rods[i][5]) / 2.0;
    li[2] = rods[i][5] - rods[i][2];

    //  length vector of the rod
    cross(B_rate, ri, Vi);
    Vi[0] = (u + Vi[0]) + VProp_Axial[18 + i % 32];
    Vi[1] += v;
    Vi[2] += w;
    cross(Vi, li, V_l_l);
    Mz_ts[0] = V_l_l[0];
    Mz_ts[1] = V_l_l[1];
    Mz_ts[2] = V_l_l[2];
    vthr_y = b_norm(Vi);
    if (!(vthr_y > 1.0E-6)) {
      vthr_y = 1.0E-6;
    }

    a = b_norm(Mz_ts) / vthr_y / b_norm(li);
    cross(Mz_ts, li, V_l_l);
    if ((V_l_l[0] * Vi[0] + V_l_l[1] * Vi[1]) + V_l_l[2] * Vi[2] < 0.0) {
      //  Vi pointing x forward, Nui pointing backward
      psiT = b_norm(li);
      vthr_y = b_norm(Vi) * (psiT * psiT);
      if (!(vthr_y > 1.0E-6)) {
        vthr_y = 1.0E-6;
      }

      V_l_l[0] /= vthr_y;
      V_l_l[1] /= vthr_y;
      V_l_l[2] /= vthr_y;
    } else {
      vthr_y = b_norm(li);
      vthr_y = b_norm(Vi) * (vthr_y * vthr_y);
      if (!(vthr_y > 1.0E-6)) {
        vthr_y = 1.0E-6;
      }

      V_l_l[0] = -V_l_l[0] / vthr_y;
      V_l_l[1] = -V_l_l[1] / vthr_y;
      V_l_l[2] = -V_l_l[2] / vthr_y;
    }

    psiT = b_norm(Vi);
    a = 0.6125 * (psiT * psiT) * b_norm(li) * rods[i][6] * 1.1 * (a * a);
    d0 = a * V_l_l[0];
    V_l_l[0] = d0;
    D_rods_idx_0 += d0;
    d0 = a * V_l_l[1];
    V_l_l[1] = d0;
    D_rods_idx_1 += d0;
    d0 = a * V_l_l[2];
    V_l_l[2] = d0;
    D_rods_idx_2 += d0;
    cross(ri, V_l_l, Vi);
    M_rods_idx_0 += Vi[0];
    M_rods_idx_1 += Vi[1];
    M_rods_idx_2 += Vi[2];
  }

  // --------------------------------------------------------------------------
  //                    Aerodynamic Forces & Moments
  // --------------------------------------------------------------------------
  //  x-direction forces
  power(vws_xcomp, dv16);
  power(vws_zcomp, Mz_ws);
  for (i = 0; i < 7; i++) {
    vws_xcomp[i] = cosf(a_ws[i]);
    a_ws[i] = sinf(a_ws[i]);
    dv16[i] += Mz_ws[i];
  }

  c_sqrt(dv16);
  power(dv16, My_ws);
  power(vwp_xcomp, dv16);
  power(vwp_zcomp, Mz_ws);
  for (i = 0; i < 7; i++) {
    d0 = Fx_ws_tmp[i] * My_ws[i];
    My_ws[i] = d0;
    Fx_ws[i] = d0 * (CL_ws[i] * a_ws[i] - CD_ws[i] * vws_xcomp[i]);
    Mx_wp[i] = cosf(a_wp[i]);
    a_wp[i] = sinf(a_wp[i]);
    dv16[i] += Mz_ws[i];
  }

  c_sqrt(dv16);
  power(dv16, vwp_zcomp);
  for (i = 0; i < 7; i++) {
    d0 = Fx_ws_tmp[i] * vwp_zcomp[i];
    vwp_zcomp[i] = d0;
    vwp_xcomp[i] = d0 * (CL_wp[i] * a_wp[i] - CD_wp[i] * Mx_wp[i]);
  }

  b_power(vts_xcomp, Mz_ts);
  li[0] = cosf(a_ts[0]);
  a_ts[0] = sinf(a_ts[0]);
  Vi[0] = Mz_ts[0];
  li[1] = cosf(a_ts[1]);
  a_ts[1] = sinf(a_ts[1]);
  Vi[1] = Mz_ts[1];
  li[2] = cosf(a_ts[2]);
  a_ts[2] = sinf(a_ts[2]);
  Vi[2] = Mz_ts[2];
  b_power(vts_zcomp, Mz_ts);
  Vi[0] += Mz_ts[0];
  Vi[1] += Mz_ts[1];
  Vi[2] += Mz_ts[2];
  d_sqrt(Vi);
  b_power(Vi, vts_xcomp);
  b_power(vtp_xcomp, Mz_ts);
  d0 = 0.0047858517437500006 * vts_xcomp[0];
  vts_xcomp[0] = d0;
  B_rate[0] = d0 * (CL_ts[0] * a_ts[0] - vthr_x * li[0]);
  ri[0] = cosf(a_tp[0]);
  a_tp[0] = sinf(a_tp[0]);
  Vi[0] = Mz_ts[0];
  d0 = 0.004557457905 * vts_xcomp[1];
  vts_xcomp[1] = d0;
  B_rate[1] = d0 * (CL_ts[1] * a_ts[1] - CD_ts_idx_1 * li[1]);
  ri[1] = cosf(a_tp[1]);
  a_tp[1] = sinf(a_tp[1]);
  Vi[1] = Mz_ts[1];
  d0 = 0.0054382547100000015 * vts_xcomp[2];
  B_rate[2] = d0 * (d1 * a_ts[2] - d2 * li[2]);
  ri[2] = cosf(a_tp[2]);
  a_tp[2] = sinf(a_tp[2]);
  Vi[2] = Mz_ts[2];
  b_power(vtp_zcomp, Mz_ts);
  Vi[0] += Mz_ts[0];
  Vi[1] += Mz_ts[1];
  Vi[2] += Mz_ts[2];
  d_sqrt(Vi);
  b_power(Vi, vtp_zcomp);
  a = 0.0047858517437500006 * vtp_zcomp[0];
  vtp_zcomp[0] = a;
  V_l_l[0] = a * (CL_tp[0] * a_tp[0] - CD_tp_idx_0 * ri[0]);
  a = 0.004557457905 * vtp_zcomp[1];
  vtp_zcomp[1] = a;
  V_l_l[1] = a * (CL_tp[1] * a_tp[1] - CD_tp_idx_1 * ri[1]);
  a = 0.0054382547100000015 * vtp_zcomp[2];
  V_l_l[2] = a * (CFy * a_tp[2] - vthr_x_tmp * ri[2]);
  c_power(vr_xcomp, dv17);
  c_power(vr_ycomp, CtrlDef);
  vthr_y = cosf(a_r[0]);
  a_r[0] = sinf(a_r[0]);
  dv17[0] += CtrlDef[0];
  vthr_z = cosf(a_r[1]);
  a_r[1] = sinf(a_r[1]);
  dv17[1] += CtrlDef[1];
  psiT = cosf(a_r[2]);
  a_r[2] = sinf(a_r[2]);
  dv17[2] += CtrlDef[2];
  CFx = cosf(a_r[3]);
  a_r[3] = sinf(a_r[3]);
  dv17[3] += CtrlDef[3];
  e_sqrt(dv17);
  c_power(dv17, vr_xcomp);
  c_power(vB_xcomp, dv17);
  c_power(vB_ycomp, CtrlDef);
  d11 = 0.0056893324250000014 * vr_xcomp[0];
  vr_xcomp[0] = d11;
  vr_ycomp[0] = d11 * (CL_r[0] * a_r[0] - CD_r_idx_0 * vthr_y);
  Mx_r[0] = cosf(a_B[0]);
  a_B[0] = sinf(a_B[0]);
  dv17[0] += CtrlDef[0];
  d11 = 0.00734973080625 * vr_xcomp[1];
  vr_xcomp[1] = d11;
  vr_ycomp[1] = d11 * (CL_r[1] * a_r[1] - CD_r_idx_1 * vthr_z);
  Mx_r[1] = cosf(a_B[1]);
  a_B[1] = sinf(a_B[1]);
  dv17[1] += CtrlDef[1];
  d11 = 0.007195253038750001 * vr_xcomp[2];
  vr_xcomp[2] = d11;
  vr_ycomp[2] = d11 * (CL_r[2] * a_r[2] - CD_r_idx_2 * psiT);
  Mx_r[2] = cosf(a_B[2]);
  a_B[2] = sinf(a_B[2]);
  dv17[2] += CtrlDef[2];
  d11 = 0.003293678325 * vr_xcomp[3];
  vr_ycomp[3] = d11 * (d5 * a_r[3] - d6 * CFx);
  Mx_r[3] = cosf(a_B[3]);
  a_B[3] = sinf(a_B[3]);
  dv17[3] += CtrlDef[3];
  e_sqrt(dv17);
  c_power(dv17, vB_xcomp);

  //  y-direction forces
  //  INCLUDE FRICTION DRAG
  //  INCLUDE FRICTION DRAG
  //  INCLUDE FRICTION DRAG
  //  INCLUDE FRICTION DRAG
  d12 = 0.01409860934125 * vB_xcomp[0];
  vB_xcomp[0] = d12;
  vB_ycomp[0] = d12 * (CL_B[0] * a_B[0] - CD_B_idx_0 * Mx_r[0]);
  CL_r[0] = vr_xcomp[0] * (-CL_r[0] * vthr_y - CD_r_idx_0 * a_r[0]);
  CL_B[0] = d12 * (-CL_B[0] * Mx_r[0] - CD_B_idx_0 * a_B[0]);
  d12 = 0.01409860934125 * vB_xcomp[1];
  vB_xcomp[1] = d12;
  vB_ycomp[1] = d12 * (CL_B[1] * a_B[1] - CD_B_idx_1 * Mx_r[1]);
  CL_r[1] = vr_xcomp[1] * (-CL_r[1] * vthr_z - CD_r_idx_1 * a_r[1]);
  CL_B[1] = d12 * (-CL_B[1] * Mx_r[1] - CD_B_idx_1 * a_B[1]);
  d12 = 0.014644945192500002 * vB_xcomp[2];
  vB_xcomp[2] = d12;
  vB_ycomp[2] = d12 * (CL_B[2] * a_B[2] - CD_B_idx_2 * Mx_r[2]);
  CL_r[2] = vr_xcomp[2] * (-CL_r[2] * psiT - CD_r_idx_2 * a_r[2]);
  CL_B[2] = d12 * (-CL_B[2] * Mx_r[2] - CD_B_idx_2 * a_B[2]);
  d12 = 0.014644945192500002 * vB_xcomp[3];
  vB_ycomp[3] = d12 * (d8 * a_B[3] - d9 * Mx_r[3]);
  CL_r[3] = d11 * (-d5 * CFx - d6 * a_r[3]);
  CL_B[3] = d12 * (-d8 * Mx_r[3] - d9 * a_B[3]);

  //  z-direction forces
  for (i = 0; i < 7; i++) {
    CL_ws[i] = My_ws[i] * (-CL_ws[i] * vws_xcomp[i] - CD_ws[i] * a_ws[i]);
    CL_wp[i] = vwp_zcomp[i] * (-CL_wp[i] * Mx_wp[i] - CD_wp[i] * a_wp[i]);
  }

  CL_ts[0] = vts_xcomp[0] * (-CL_ts[0] * li[0] - vthr_x * a_ts[0]);
  CL_tp[0] = vtp_zcomp[0] * (-CL_tp[0] * ri[0] - CD_tp_idx_0 * a_tp[0]);
  CL_ts[1] = vts_xcomp[1] * (-CL_ts[1] * li[1] - CD_ts_idx_1 * a_ts[1]);
  CL_tp[1] = vtp_zcomp[1] * (-CL_tp[1] * ri[1] - CD_tp_idx_1 * a_tp[1]);
  CL_ts[2] = d0 * (-d1 * li[2] - d2 * a_ts[2]);
  CL_tp[2] = a * (-CFy * ri[2] - vthr_x_tmp * a_tp[2]);

  //  INCLUDE FRICTION DRAG
  //  INCLUDE FRICTION DRAG
  //  x-direction moments
  for (i = 0; i < 7; i++) {
    vws_zcomp[i] = b[i] * CL_ws[i];
    Mx_wp[i] = c_b[i] * CL_wp[i];
  }

  //  y-direction moments
  for (i = 0; i < 7; i++) {
    My_ws[i] = (-0.00661 * Fx_ws[i] - b_b[i] * CL_ws[i]) + My_ws[i] * Cw[i] *
      CM_ws[i];
    vwp_zcomp[i] = (-0.00661 * vwp_xcomp[i] - b_b[i] * CL_wp[i]) + vwp_zcomp[i] *
      Cw[i] * CM_wp[i];
  }

  vts_xcomp[0] = (-0.00661 * B_rate[0] - -0.47126999999999997 * CL_ts[0]) +
    vts_xcomp[0] * 0.11741 * deltaT;
  vtp_zcomp[0] = (-0.00661 * V_l_l[0] - -0.47126999999999997 * CL_tp[0]) +
    vtp_zcomp[0] * 0.11741 * CM_tp_idx_0;
  vts_xcomp[1] = (-0.00661 * B_rate[1] - -0.47972000000000004 * CL_ts[1]) +
    vts_xcomp[1] * 0.14564 * CM_ts_idx_1;
  vtp_zcomp[1] = (-0.00661 * V_l_l[1] - -0.47972000000000004 * CL_tp[1]) +
    vtp_zcomp[1] * 0.14564 * CM_tp_idx_1;

  //  z-direction moments
  for (i = 0; i < 7; i++) {
    Mz_ws[i] = 0.0 - b[i] * Fx_ws[i];
    vws_xcomp[i] = 0.0 - c_b[i] * vwp_xcomp[i];
  }

  //  WHY NEGATIVE SIGN???
  vr_xcomp[0] = (-0.49394000000000005 * CL_r[0] - 0.0 * vr_ycomp[0]) - vr_xcomp
    [0] * 0.20062000000000002 * CM_r_idx_0;
  vB_xcomp[0] = (0.03491 * CL_B[0] - 0.0 * vB_ycomp[0]) - vB_xcomp[0] * 0.64171 *
    CM_B_idx_0;
  vr_xcomp[1] = (-0.49581 * CL_r[1] - 0.0 * vr_ycomp[1]) - vr_xcomp[1] *
    0.19294999999999998 * CM_r_idx_1;
  vB_xcomp[1] = (0.03491 * CL_B[1] - 0.0 * vB_ycomp[1]) - vB_xcomp[1] * 0.64171 *
    CM_B_idx_1;
  vr_xcomp[2] = (-0.50414 * CL_r[2] - 0.0 * vr_ycomp[2]) - vr_xcomp[2] * 0.16949
    * CM_r_idx_2;
  vB_xcomp[2] = (0.03491 * CL_B[2] - 0.0 * vB_ycomp[2]) - vB_xcomp[2] * 0.64171 *
    CM_B_idx_2;

  //  WHY NEGATIVE SIGN???
  //  Total forces and moments in x, y and z directions
  vthr_z = Fx_ws[0];
  psiT = vwp_xcomp[0];
  for (i = 0; i < 6; i++) {
    vthr_z += Fx_ws[i + 1];
    psiT += vwp_xcomp[i + 1];
  }

  *Fx = ((((((vthr_z + psiT) + ((B_rate[0] + B_rate[1]) + B_rate[2])) + ((V_l_l
              [0] + V_l_l[1]) + V_l_l[2])) + (((vr_ycomp[0] + vr_ycomp[1]) +
             vr_ycomp[2]) + vr_ycomp[3])) + (((vB_ycomp[0] + vB_ycomp[1]) +
            vB_ycomp[2]) + vB_ycomp[3])) + ThrForce_idx_0) + D_rods_idx_0;
  *Fy = (((((CL_r[0] + CL_r[1]) + CL_r[2]) + CL_r[3]) + (((CL_B[0] + CL_B[1]) +
            CL_B[2]) + CL_B[3])) + ThrForce_idx_1) + D_rods_idx_1;
  vthr_z = CL_ws[0];
  psiT = CL_wp[0];
  for (i = 0; i < 6; i++) {
    vthr_z += CL_ws[i + 1];
    psiT += CL_wp[i + 1];
  }

  *Fz = ((((vthr_z + psiT) + ((CL_ts[0] + CL_ts[1]) + CL_ts[2])) + ((CL_tp[0] +
            CL_tp[1]) + CL_tp[2])) + ThrForce_idx_2) + D_rods_idx_2;
  vthr_z = vws_zcomp[0];
  psiT = Mx_wp[0];
  for (i = 0; i < 6; i++) {
    vthr_z += vws_zcomp[i + 1];
    psiT += Mx_wp[i + 1];
  }

  *Mx = (((((((vthr_z + psiT) + ((0.03953 * CL_ts[0] + 0.09333 * CL_ts[1]) +
    0.15156999999999998 * CL_ts[2])) + ((-0.03953 * CL_tp[0] + -0.09333 * CL_tp
    [1]) + -0.15156999999999998 * CL_tp[2])) + ((((0.0 - 0.016929999999999997 *
    CL_r[0]) + (0.0 - -0.03657 * CL_r[1])) + (0.0 - -0.10124 * CL_r[2])) + (0.0
              - -0.15517 * CL_r[3]))) + ((((0.0 - 0.047200000000000006 * CL_B[0])
              + (0.0 - 0.01133 * CL_B[1])) + (0.0 - -0.025240000000000002 *
              CL_B[2])) + (0.0 - -0.0625 * CL_B[3]))) + (0.0 * ThrForce_idx_2 -
           -0.00661 * ThrForce_idx_1)) + MX) + M_rods_idx_0;
  vthr_z = My_ws[0];
  psiT = vwp_zcomp[0];
  for (i = 0; i < 6; i++) {
    vthr_z += My_ws[i + 1];
    psiT += vwp_zcomp[i + 1];
  }

  vthr_y = wOut * 2.0 * M_PI / 60.0;
  *My = ((((((((vthr_z + psiT) + ((vts_xcomp[0] + vts_xcomp[1]) + ((-0.00661 *
    B_rate[2] - -0.47514999999999996 * CL_ts[2]) + d0 * 0.14316 * d3))) +
              ((vtp_zcomp[0] + vtp_zcomp[1]) + ((-0.00661 * V_l_l[2] -
    -0.47514999999999996 * CL_tp[2]) + a * 0.14316 * d4))) +
             (((0.016929999999999997 * vr_ycomp[0] + -0.03657 * vr_ycomp[1]) +
               -0.10124 * vr_ycomp[2]) + -0.15517 * vr_ycomp[3])) +
            (((0.047200000000000006 * vB_ycomp[0] + 0.01133 * vB_ycomp[1]) +
              -0.025240000000000002 * vB_ycomp[2]) + -0.0625 * vB_ycomp[3])) + (
            -0.00661 * ThrForce_idx_0 - 0.24 * ThrForce_idx_2)) + 0.5 * (1.225 *
           (b_a * b_a) * 0.001057227821024 * (-CMy * b_ThrForce_idx_1_tmp + CMz *
            ThrForce_idx_1_tmp))) + -5.0E-5 * vthr_y * r) + M_rods_idx_1;
  vthr_z = Mz_ws[0];
  psiT = vws_xcomp[0];
  for (i = 0; i < 6; i++) {
    vthr_z += Mz_ws[i + 1];
    psiT += vws_xcomp[i + 1];
  }

  *Mz = ((((((((vthr_z + psiT) + (((-0.0 - 0.03953 * B_rate[0]) + (-0.0 -
    0.09333 * B_rate[1])) + (-0.0 - 0.15156999999999998 * B_rate[2]))) + (((-0.0
    - -0.03953 * V_l_l[0]) + (-0.0 - -0.09333 * V_l_l[1])) + (-0.0 -
    -0.15156999999999998 * V_l_l[2]))) + (((vr_xcomp[0] + vr_xcomp[1]) +
    vr_xcomp[2]) + ((-0.52076 * CL_r[3] - 0.0 * vr_ycomp[3]) - d11 * 0.14114 *
                    d7))) + (((vB_xcomp[0] + vB_xcomp[1]) + vB_xcomp[2]) +
             ((0.03491 * CL_B[3] - 0.0 * vB_ycomp[3]) - d12 * 0.64171 * d10))) +
           (0.24 * ThrForce_idx_1 - 0.0 * ThrForce_idx_0)) + 0.5 * (1.225 * (c_a
            * c_a) * 0.001057227821024 * (-CMy * ThrForce_idx_1_tmp - CMz *
            b_ThrForce_idx_1_tmp))) + 5.0E-5 * vthr_y * q) + M_rods_idx_2;
}

//
// Arguments    : void
// Return Type  : void
//

//
// File trailer for McFoamy_FM_v2.cpp
//
// [EOF]
//
