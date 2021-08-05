//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: McFoamy_FM_v3.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 05-Aug-2021 13:41:18
//

// Include Files
// #include "McFoamy_FM_v3.h"
#include "flappedAirfoil_LAR_v1_FlatPlate.h"
#include "interp2.h"
#include "rt_nonfinite.h"
#include "sign.h"
#include "rt_defines.h"
#include <cmath>


// real_T rtNaN{std::numeric_limits<real_T>::quiet_NaN()};
// real_T rtInf{std::numeric_limits<real_T>::infinity()};
// real_T rtMinusInf{-std::numeric_limits<real_T>::infinity()};
// real32_T rtNaNF{std::numeric_limits<real32_T>::quiet_NaN()};
// real32_T rtInfF{std::numeric_limits<real32_T>::infinity()};
// real32_T rtMinusInfF{-std::numeric_limits<real32_T>::infinity()};

#include "allMethods.hpp"

// Function Declarations
static double rt_atan2d_snf(double u0, double u1);

// Function Definitions
//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  if (std::isnan(u0) || std::isnan(u1)) {
    y = rtNaN;
  } else if (std::isinf(u0) && std::isinf(u1)) {
    int b_u0;
    int b_u1;
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
    y = std::atan2(static_cast<double>(b_u0), static_cast<double>(b_u1));
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = std::atan2(u0, u1);
  }
  return y;
}

//
// Actuator limits
//
// Arguments    : double Ail_def
//                double Elev_def
//                double Rud_def
//                double Thr_com
//                double v_u
//                double v_v
//                double v_w
//                double w_p
//                double w_q
//                double w_r
//                double *Fx
//                double *Fy
//                double *Fz
//                double *Mx
//                double *My
//                double *Mz
// Return Type  : void
//
void McFoamy_FM_v3(double Ail_def, double Elev_def, double Rud_def,
                   double Thr_com, double v_u, double v_v, double v_w,
                   double w_p, double w_q, double w_r,  double *Fx, double *Fy, double *Fz, double *Mx, double *My,
                   double *Mz)
{
  Ail_def *= -57.2958;
  Elev_def *= 57.2958;
  Rud_def *= 57.2958;

  static const double Cts_Geom[154]{
      7.0,   260.41,   152.4,  431.8,  0.0,     0.0,    0.0,
      3.0,   152.73,   139.7,  181.61, 0.0,     0.0,    0.0,
      4.0,   204.01,   128.55, 215.9,  0.0,     0.0,    0.0,
      4.0,   641.71,   641.71, 146.27, 0.0,     0.0,    0.0,
      57.26, 14499.66, 253.4,  0.0,    -217.07, 33.38,  0.0,
      58.72, 13977.71, 237.58, 102.04, -211.39, 91.76,  0.0,
      59.29, 13198.25, 222.69, 99.23,  -215.43, 151.14, 0.0,
      62.58, 12974.4,  207.42, 96.36,  -212.18, 212.01, 0.0,
      62.58, 11992.21, 191.74, 93.41,  -208.85, 274.55, 0.0,
      63.25, 11122.51, 175.97, 90.44,  -205.5,  337.42, 0.0,
      61.9,  9913.9,   160.28, 87.49,  -202.26, 400.21, 0.0,
      66.55, 7241.76,  117.41, 81.87,  -711.27, 39.53,  0.0,
      51.09, 7452.58,  145.64, 111.49, -719.72, 93.33,  0.0,
      62.02, 8892.67,  143.16, 143.16, -715.15, 151.57, 0.0,
      46.3,  9291.95,  200.62, 120.36, -733.94, 0.0,    23.54,
      62.19, 11997.38, 192.95, 112.69, -735.81, 0.0,    -29.96,
      69.31, 11715.82, 169.49, 103.43, -744.14, 0.0,    -94.63,
      38.1,  5298.95,  141.14, 141.14, -760.76, 0.0,    -148.56,
      35.87, 23018.14, 641.71, 0.0,    -205.09, 0.0,    53.81,
      35.87, 23018.14, 641.71, 0.0,    -205.09, 0.0,    17.94,
      37.26, 23910.11, 641.71, 0.0,    -205.09, 0.0,    -18.63,
      37.26, 23910.11, 641.71, 0.0,    -205.09, 0.0,    -55.89};
  static const double b_dv[110]{
      0.156336454,  0.156336454,  0.156336454,  0.156336454,  0.156336454,
      0.156336454,  0.156336454,  0.156336454,  0.156336454,  0.156336454,
      0.142939755,  0.143274298,  0.144236696,  0.145715204,  0.147592844,
      0.149730239,  0.151960378,  0.154120349,  0.156099101,  0.157836404,
      0.12449467,   0.125294853,  0.127662714,  0.131491476,  0.136533063,
      0.142227959,  0.148100411,  0.153647038,  0.158176177,  0.161723393,
      0.103654093,  0.104978752,  0.10890636,   0.115283154,  0.123846679,
      0.133924798,  0.144287776,  0.153875413,  0.161071479,  0.163852071,
      0.080437634,  0.082320505,  0.087900516,  0.096990697,  0.109258684,
      0.124116893,  0.139693524,  0.153953862,  0.163349211,  0.165797264,
      0.054966543,  0.057426888,  0.064729471,  0.076657686,  0.092816061,
      0.112598271,  0.133992265,  0.153533548,  0.165846339,  0.168474477,
      0.027444178,  0.03049664,   0.039588216,  0.054473185,  0.074693054,
      0.099537509,  0.127172755,  0.152570626,  0.168638349,  0.171699702,
      -0.001874328, 0.001785855,  0.012726044,  0.030695378,  0.05514811,
      0.085260174,  0.119365045,  0.151139537,  0.171659606,  0.175259993,
      -0.032728784, -0.02842523,  -0.015539759, 0.005645211,  0.034482474,
      0.070037437,  0.11072865,   0.149362748,  0.174874022,  0.179027738,
      -0.064814458, -0.059769732, -0.044818343, -0.020318541, 0.013005062,
      0.054122603,  0.101419857,  0.147360112,  0.17827724,   0.182929631,
      -0.097502535, -0.091725529, -0.074619748, -0.046776089, -0.008957689,
      0.037754768,  0.091641937,  0.145232828,  0.181871525,  0.18696977};
  static const double b_dv1[110]{
      0.0,         0.0,         0.0,         0.0,         0.0,
      0.0,         0.0,         0.0,         0.0,         0.0,
      0.0,         0.000304035, 0.000593774, 0.000855707, 0.001081293,
      0.001265997, 0.001410344, 0.001515466, 0.001581189, 0.00160638,
      0.0,         0.000672049, 0.001308186, 0.001875984, 0.002343891,
      0.002685232, 0.002915524, 0.003073033, 0.003203695, 0.003296003,
      0.0,         0.001105525, 0.002139399, 0.003046099, 0.003777072,
      0.004271303, 0.004529052, 0.004674445, 0.004870829, 0.005159087,
      0.0,         0.001600104, 0.003086457, 0.004372415, 0.005383801,
      0.006042688, 0.006286839, 0.006333996, 0.006596278, 0.007149694,
      0.0,         0.00214306,  0.004124254, 0.005829393, 0.007152135,
      0.007991019, 0.008217389, 0.008079318, 0.008340642, 0.00921631,
      0.0,         0.002695169, 0.005206331, 0.007370573, 0.009046337,
      0.010092937, 0.010327588, 0.009939438, 0.010095901, 0.01134185,
      0.0,         0.003243617, 0.0062742,   0.008920992, 0.010988567,
      0.012301024, 0.012588497, 0.01193205,  0.011869291, 0.013520997,
      0.0,         0.003737499, 0.007197636, 0.010310416, 0.012837704,
      0.014518175, 0.014938663, 0.014050094, 0.013673189, 0.015752967,
      0.0,         0.004001698, 0.007738064, 0.011301405, 0.014405295,
      0.016623841, 0.017317246, 0.016267158, 0.015518142, 0.018038805,
      0.0,         0.003755281, 0.007554928, 0.011576152, 0.015455208,
      0.018472368, 0.019651635, 0.018549771, 0.017413651, 0.020375896};
  static const double b_dv2[110]{
      -0.009379255, -0.009379255, -0.009379255, -0.009379255, -0.009379255,
      -0.009379255, -0.009379255, -0.009379255, -0.009379255, -0.009379255,
      -0.009414602, -0.009414258, -0.009412543, -0.009408081, -0.009399845,
      -0.009386461, -0.009367587, -0.009343565, -0.009317827, -0.009292775,
      -0.00922517,  -0.009234779, -0.009259831, -0.009292432, -0.0093192,
      -0.009327092, -0.009306845, -0.009256742, -0.009186392, -0.009135945,
      -0.008715215, -0.008749876, -0.008846651, -0.008982547, -0.009121875,
      -0.009224484, -0.009256742, -0.009206982, -0.009150015, -0.009359351,
      -0.007778011, -0.007856254, -0.008075885, -0.008395035, -0.008745071,
      -0.009041916, -0.009209727, -0.009213159, -0.00927596,  -0.009773903,
      -0.006305113, -0.0064465,   -0.006846983, -0.007438956, -0.008111231,
      -0.008720706, -0.009127023, -0.00925331,  -0.009427299, -0.01020081,
      -0.004198032, -0.004423154, -0.005068319, -0.006031948, -0.007148975,
      -0.008203888, -0.008967104, -0.009300668, -0.009571088, -0.010605068,
      -0.001374064, -0.001706255, -0.002664393, -0.004110867, -0.005807856,
      -0.007449251, -0.00869943,  -0.009331211, -0.009705269, -0.010994226,
      0.002229593,  0.001762192,  0.00040563,   -0.001648602, -0.004067627,
      -0.006438264, -0.008307869, -0.009330867, -0.009834645, -0.011382011,
      0.006644511,  0.00599214,   0.004131457,  0.001338374,  -0.001941328,
      -0.005174359, -0.007785904, -0.009294491, -0.009963678, -0.011781807,
      0.011790387,  0.01091461,   0.008431755,  0.004772847,  0.000516475,
      -0.003684646, -0.007140396, -0.009222768, -0.010097172, -0.012199449};
  static const double b_dv3[110]{
      0.0,         0.0,         0.0,         0.0,         0.0,
      0.0,         0.0,         0.0,         0.0,         0.0,
      0.0,         0.000624918, 0.001228902, 0.001786557, 0.002280383,
      0.002695622, 0.003025754, 0.003265632, 0.003407706, 0.003447857,
      0.0,         0.00118875,  0.002355195, 0.003475311, 0.004508604,
      0.005373056, 0.006031948, 0.006465718, 0.006690153, 0.006735452,
      0.0,         0.001686008, 0.003352112, 0.004976005, 0.006530234,
      0.007916309, 0.008970193, 0.009615701, 0.009827095, 0.009570745,
      0.0,         0.00210811,  0.004209014, 0.006291729, 0.008331891,
      0.010270131, 0.01181784,  0.012740975, 0.012792794, 0.012191899,
      0.0,         0.002447508, 0.004910802, 0.007398462, 0.009900534,
      0.012377212, 0.01452513,  0.015837766, 0.015769475, 0.014745447,
      0.0,         0.002695622, 0.005441004, 0.008275611, 0.011210425,
      0.014212843, 0.017039214, 0.018875875, 0.018797288, 0.017259188,
      0.0,         0.002844902, 0.005784863, 0.008902245, 0.012237884,
      0.015771877, 0.019318224, 0.021818238, 0.021878979, 0.019738611,
      0.0,         0.002886426, 0.005919044, 0.009246104, 0.012952713,
      0.017029948, 0.021326814, 0.024630882, 0.025002195, 0.022188521,
      0.0,         0.00279274,  0.005797904, 0.009258801, 0.013314074,
      0.017960633, 0.023040962, 0.027285666, 0.028151834, 0.024615439,
      0.0,         0.002507907, 0.005340455, 0.008857975, 0.013262255,
      0.018530643, 0.024451059, 0.029763374, 0.031305935, 0.027038239};
  static const double b_dv4[110]{
      0.0,         0.0,         0.0,          0.0,          0.0,
      0.0,         0.0,         0.0,          0.0,          0.0,
      0.0,         0.000466715, 0.000928969,  0.001379898,  0.001812295,
      0.002217926, 0.002588209, 0.002912507,  0.003173662,  0.003355543,
      0.0,         0.000650999, 0.001316754,  0.002009277,  0.002733028,
      0.003466731, 0.004186021, 0.004849031,  0.005419041,  0.005744369,
      0.0,         0.000651342, 0.001337687,  0.002091295,  0.002939961,
      0.003887461, 0.004884035, 0.005825014,  0.00627663,   0.005094056,
      0.0,         0.000542556, 0.001141049,  0.001847642,  0.002706947,
      0.003752251, 0.004944433, 0.006105387,  0.005895365,  0.003495558,
      0.0,         0.000378176, 0.000833567,  0.001435835,  0.002244007,
      0.003306126, 0.004618762, 0.005966058,  0.005365163,  0.001968096,
      0.0,         0.000193549, 0.00048456,   0.000959511,  0.001687037,
      0.002723076, 0.004095081, 0.005597491,  0.004922813,  0.000671589,
      0.0,         7.5498E-6,   0.000133151,  0.000477353,  0.001115998,
      0.002106738, 0.003494871, 0.005116363,  0.004597142,  -0.000392246,
      0.0,         -0.00017193, -0.000201786, 1.9904E-5,    0.000571383,
      0.001510303, 0.002888485, 0.004593367,  0.004365844,  -0.001249149,
      0.0,         -0.00033528, -0.000507209, -0.000396021, 7.41253E-5,
      0.000958825, 0.002311955, 0.004069343,  0.004200435,  -0.001923827,
      0.0,         -0.00047049, -0.000764589, -0.000754637, -0.000362734,
      0.000464313, 0.001782782, 0.003567967,  0.004069686,  -0.002446135};
  double VProp_Axial[18];
  double CL_wp[7];
  double CL_ws[7];
  double CM_wp[7];
  double CM_ws[7];
  double Cw[7];
  double Cwf[7];
  double a_wp[7];
  double a_ws[7];
  double vwp_xcomp[7];
  double vwp_xz[7];
  double vwp_zcomp[7];
  double vws_xcomp[7];
  double vws_xz[7];
  double vws_zcomp[7];
  double xw[7];
  double ywp[7];
  double yws[7];
  double zw[7];
  double vB_ycomp[4];
  double vr_ycomp[4];
  double AilDefOut;
  double CB_idx_0;
  double CB_idx_1;
  double CB_idx_2;
  double CB_idx_3;
  double CFY_tmp;
  double CFx;
  double CFy;
  double CL_B_idx_0;
  double CL_B_idx_1;
  double CL_B_idx_2;
  double CL_r_idx_0;
  double CL_r_idx_1;
  double CL_r_idx_2;
  double CL_tp_idx_0;
  double CL_tp_idx_1;
  double CL_ts_idx_0;
  double CL_ts_idx_1;
  double CM_B_idx_0;
  double CM_B_idx_1;
  double CM_B_idx_2;
  double CM_r_idx_0;
  double CM_r_idx_1;
  double CM_r_idx_2;
  double CM_tp_idx_0;
  double CM_tp_idx_1;
  double CM_ts_idx_0;
  double CM_ts_idx_1;
  double CMx;
  double CMy;
  double CMz;
  double Cr_idx_0;
  double Cr_idx_1;
  double Cr_idx_2;
  double Cr_idx_3;
  double Ct_idx_0;
  double Ct_idx_1;
  double Ct_idx_2;
  double Ctf_idx_0;
  double Ctf_idx_1;
  double Ctf_idx_2;
  double CtrlDef_idx_3;
  double ElevDefOut;
  double F_Thr_idx_0;
  double F_Thr_idx_1;
  double F_Thr_idx_2;
  double Fx_ts_tmp_idx_0;
  double Fx_ts_tmp_idx_1;
  double Fx_ts_tmp_idx_2;
  double MX_thr;
  double MX_thr_tmp;
  double RAilDef;
  double RudDefOut;
  double ThrRPMOut;
  double W_Thr_SI;
  double a_B_idx_0;
  double a_r_idx_0;
  double a_tp_idx_0;
  double a_tp_idx_1;
  double a_tp_idx_2;
  double a_ts_idx_0;
  double a_ts_idx_1;
  double a_ts_idx_2;
  double b_CFY_tmp;
  double b_y_idx_0;
  double d;
  double d1;
  double d10;
  double d11;
  double d12;
  double d13;
  double d14;
  double d15;
  double d16;
  double d2;
  double d3;
  double d4;
  double d5;
  double d6;
  double d7;
  double d8;
  double d9;
  double psiT;
  double vB_xy_idx_0;
  double vB_xy_idx_1;
  double vB_xy_idx_2;
  double vB_xy_idx_3;
  double vr_xy_idx_0;
  double vr_xy_idx_1;
  double vr_xy_idx_2;
  double vr_xy_idx_3;
  double vthr_x;
  double vthr_z;
  double vtp_xcomp_idx_0;
  double vtp_xcomp_idx_1;
  double vtp_xz_idx_0;
  double vtp_xz_idx_1;
  double vtp_xz_idx_2;
  double vtp_zcomp_idx_0;
  double vtp_zcomp_idx_1;
  double vts_xcomp_idx_0;
  double vts_xcomp_idx_1;
  double vts_xcomp_idx_2;
  double vts_xz_idx_1;
  double vts_xz_idx_2;
  double vts_zcomp_idx_0;
  double vts_zcomp_idx_1;
  double xB_idx_0;
  double xB_idx_1;
  double xB_idx_2;
  double xB_idx_3;
  double xr_idx_0;
  double xr_idx_1;
  double xr_idx_2;
  double xr_idx_3;
  double xt_idx_0;
  double xt_idx_1;
  double xt_idx_2;
  double yB_idx_0;
  double yB_idx_1;
  double yB_idx_2;
  double yB_idx_3;
  double y_idx_0;
  double y_idx_1;
  double y_idx_2;
  double y_idx_3;
  double yr_idx_0;
  double yr_idx_1;
  double yr_idx_2;
  double yr_idx_3;
  double ytp_idx_0;
  double ytp_idx_1;
  double yts_idx_0;
  double yts_idx_1;
  double zB_idx_0;
  double zB_idx_1;
  double zB_idx_2;
  double zB_idx_3;
  double zr_idx_0;
  double zr_idx_1;
  double zr_idx_2;
  double zr_idx_3;
  double zt_idx_0;
  double zt_idx_1;
  double zt_idx_2;
  int Cw_tmp;
  int i;
  //  Cbi_0 = [0  0  -1;
  //           0  1  0;
  //           1  0  0;];
  //  Aircraft and propeller geometry
  // MCFOAMY GEOMETRY
  //         x     y    z
  //  CG = [-238.58, 0, 5.89]; %measured from propeller plane
  // hand measured
  // CG = [-370, 0, 5.89];
  // CG = [-0.309419,0,0.004097]*1e3;    % mm - McFoamy with VTOL gear
  //  %
  //  -------------------------------------------------------------------------
  //  % Corners for contact dynamics
  //  %
  //  -------------------------------------------------------------------------
  //  % Geometry, McFoamy corners [m] of the rigid body with respect to the CG
  //  % n=8;
  //  % R=0.127;    % radius of the propeller [m]
  //  distProp=0.04;
  //  % prop=zeros(n,3);
  //  % for i=1:n;
  //  %     prop(i,:)=R*[0,cos(i*2*pi/n), sin(i*2*pi/n)];
  //  % end
  //  c=[-.14, .408,0;
  //     -.14, -.408,0;
  //     -.272, .408, 0;
  //     -.272, -.408, 0;
  //     -.77, .178, 0;
  //     -.77, -.178, 0;
  //     -.89, 0, .06;
  //     -.77, 0, .07;
  //     -.87, 0, -.17;
  //     -.77, 0, -.17;
  //     -.12, .09, .16;
  //     -.12, -.09, .16;
  //     -.94, .2, .2;    % VTOL gears (not visible)
  //     -.94, -.2, .2;
  //     -.94, -.2, -.2;
  //     -.94, .2, -.2;];
  //  c(:,1)=c(:,1)-distProp*ones(size(c,1),1);
  //  corners=[c]; % corners=[prop;c];
  //  corners=corners-ones(size(corners,1),1)*(CG.*10^(-3));
  //  %
  //  -------------------------------------------------------------------------
  //  % Carbon fiber rods for drag calculation
  //  %
  //  -------------------------------------------------------------------------
  //
  //  % VTOL gear
  //  a1=[-562, 39.2, 0]; % attachement points [mm]
  //  a2=[-562, -39.2, 0];
  //  b1=[-715, 117.3, 0];
  //  b2=[-715, -117.3, 0];
  //  r=200;                  % half length of the rear square
  //  c1=[-972, r, -r];
  //  c2=[-972, -r, -r];
  //  c3=[-972, r, r];
  //  c4=[-972, -r, r];
  //  c5=[-972,0,0];
  //
  //  d1=2;   % small rod diameter
  //  d2=2.5; % big rod diameter
  //          %  pt1, pt2, diameter [mm]
  //  land_gear=[a1,c1,d2;    % [mm] with respect to propeller plane
  //             a1,c3,d2;
  //             a2,c2,d2;
  //             a2,c4,d2;
  //             b1,c1,d1;
  //             b1,c3,d1;
  //             b2,c2,d1;
  //             b2,c4,d1;
  //             c2,c1,d1;
  //             c1,c3,d1;
  //             c4,c2,d1;
  //             c3,c4,d1;
  //             c1,c5,d1;
  //             c2,c5,d1;
  //             c3,c5,d1;
  //             c4,c5,d1;
  //             a1,b1,d1;
  //             a2,b2,d1];
  //
  //  rod=[-80, 48, 0, -255, 0 -80;    % structural rods [mm] (7x2 top)
  //      -255, 0 -80, -328, 54, 0;
  //      -328, 54, 0, -425, 0, -85;
  //      -425, 0, -85, -540 38, 0;
  //      -540 38, 0, -645, 0, -60;
  //      -645, 0, -60, -695, 115, 0;
  //      -745, 0, -93, -720, 115, 0;
  //      -80, 48, 0, -230, 0, 75;    % 9x2 rods bottom
  //      -230, 0, 75, -185, 210, 0;
  //      -230, 0, 75, -265, 210, 0;
  //      -230, 0, 75,-330, 55, 0;
  //      -330, 55, 0, -445 0, 65;
  //      -445 0, 65, -540, 35, 0;
  //      -540, 35, 0, -655, 0, 55;
  //      -655, 0, 55, -695, 115, 0;
  //      -720, 115, 0, -745, 0, 65];
  //
  //  rod=[rod; rod*diag([1 -1 1 1 -1 1])]; % mirror it for y
  //  rod=[rod, d1*ones(size(rod,1),1)];  % add the diameter in last row
  //  %rods=[rod; land_gear]; % append the landing gears in the same table
  //  rods=[rod];
  //
  //  % rod centers [mm] w/r to prop plane
  //  rod_geom=[zeros(size(rods,1),4), (rods(:,1:3)+rods(:,4:6))/2];
  //
  //  % shift it by CG
  //  rods(:,1:6)=(rods(:,1:6)-[ones(size(rods,1),1)*CG,ones(size(rods,1),1)*CG]);
  //
  //  rods=rods*1e-3;   % convert to [m]
  //  -------------------------------------------------------------------------
  //  Elements decomposition
  //  -------------------------------------------------------------------------
  //          #segs  wing root   wing tip  wing span
  //  Wing Sections
  //             Span    Area     Chord FlapChord   x     y     z
  //  Tail Sections
  //             Span    Area     Chord FlapChord   x     y     z
  //  Rudder Sections
  //             Span    Area     Chord FlapChord   x     y     z
  //  Body Sections
  //             Span    Area    Chord  FlapChord   x     y     z
  //  WING TAIL GEOMETRY
  //  WING RUDDER GEOMETRY
  //  Propeller tables
  //  Propeller tables
  //  Negative J data from Selig's database (for APC 10x4.7) % we don't use it
  //  Maximum Actuators
  // degrees
  // degrees
  // degrees
  // RPM
  //  Inertial parameters
  //  Controller Gains (MUST CLEAN THIS)
  //  Attitude control gains
  //  Steady controller Gains
  //  Xi to roll
  //  Gains for tracking
  //  Kpp = 0.5*0.015*[150 50 500];
  //  Kpd = 0.8*0.05*[6 2 3];
  //  Kpi = 0.8*0.1*[1 0.5 3];
  //  Threshold (degrees)
  if ((Ail_def >= -52.0) && (Ail_def <= 52.0)) {
    AilDefOut = Ail_def;
  } else if (Ail_def < -52.0) {
    AilDefOut = -52.0;
  } else {
    AilDefOut = 52.0;
  }
  if ((Elev_def >= -59.0) && (Elev_def <= 59.0)) {
    ElevDefOut = Elev_def;
  } else if (Elev_def < -59.0) {
    ElevDefOut = -59.0;
  } else {
    ElevDefOut = 59.0;
  }
  if ((Rud_def >= -49.0) && (Rud_def <= 49.0)) {
    RudDefOut = Rud_def;
  } else if (Rud_def < -49.0) {
    RudDefOut = -49.0;
  } else {
    RudDefOut = 49.0;
  }
  if ((Thr_com >= 0.0) && (Thr_com <= 6710.0)) {
    ThrRPMOut = Thr_com;
  } else if (Thr_com < 0.0) {
    ThrRPMOut = 0.0;
  } else {
    ThrRPMOut = 6710.0;
  }
  CtrlDef_idx_3 = 0.017453292519943295 * RudDefOut;
  // AER_lim[0] = AilDefOut;
  // AER_lim[1] = ElevDefOut;
  // AER_lim[2] = RudDefOut;
  //  Handling of system states
  //  % dummy values for testing
  //
  //  velocity_wind = [0 0 0]'; %in inertial coordinates
  //  velocity_body = [0 0 0]';
  //  omega_body = [0 0 0]';
  //  C_bi = eye(3);
  //
  //  velocity_wind = [Velocity_wind(1); Velocity_wind(2); Velocity_wind(3)];
  //  velocity_body = [Velocity_body(1); Velocity_body(2); Velocity_body(3)];
  //  omega_body = [Omega_body(1); Omega_body(2); Omega_body(3)];
  //  C_bi = [C_bi_col(1:3)';C_bi_col(4:6)';C_bi_col(7:9)']';
  //  Airspeed calculation
  //  Vel_a = velocity_body-C_bi*velocity_wind; % airspeed
  //  wb_1 = omega_body(1); % p
  //  wb_2 = omega_body(2); % q
  //  wb_3 = omega_body(3); % r
  //  p
  //  q
  //  r
  //  Thruster model
  // Rh = 0.006; % propeller hub radius
  // Rp = 0.127; % propeller radius
  // D = 2*Rp;   % propeller dia
  vthr_x = (v_u + w_q * -0.0058899999999999994) - w_r * -0.0;
  RudDefOut = (v_v + w_r * 0.27) - w_p * -0.0058899999999999994;
  vthr_z = (v_w + w_p * -0.0) - w_q * 0.27;
  // -------------------------------------------------------------
  //  total velocity
  //  in-plane velocity
  //  azimuth angle of thruster (0 to +90 deg.)
  CMx = RudDefOut * RudDefOut;
  CFx = vthr_z * vthr_z;
  psiT = std::fmax(
      0.0, std::fmin(1.5707963267948966,
                     rt_atan2d_snf(std::sqrt(CMx + CFx), std::abs(vthr_x))));
  vthr_z = rt_atan2d_snf(vthr_z, RudDefOut);
  //  in-plane angle (0 to +/- 180)
  if (ThrRPMOut < 1716.0) {
    //  Too slow to run the motor
    ThrRPMOut = 0.0;
    CFx = 0.0;
    CFy = 0.0;
    CMx = 0.0;
    CMy = 0.0;
    CMz = 0.0;
  } else {
    //  advance ratio based on total velocity
    RudDefOut =
        std::fmax(0.0, std::fmin(1.0, std::sqrt((vthr_x * vthr_x + CMx) + CFx) /
                                          (ThrRPMOut / 60.0 * 0.254)));
    if (vthr_x >= 0.0) {
      //  Forward flight
      CFx = coder::interp2_dispatch(b_dv, RudDefOut, psiT, rtNaN);
      CFy = coder::interp2_dispatch(b_dv1, RudDefOut, psiT, rtNaN);
      CMx = coder::interp2_dispatch(b_dv2, RudDefOut, psiT, rtNaN);
      CMy = coder::interp2_dispatch(b_dv3, RudDefOut, psiT, rtNaN);
      CMz = coder::interp2_dispatch(b_dv4, RudDefOut, psiT, rtNaN);
    } else {
      //  Rearward flight
      //  A/c to Bart's paper Fig.8 & 9, for rearward flight, Cfx and Cmx have
      //  approx. static value
      CFx = 0.15633;
      CFy = 0.0;
      CMx = -0.00938;
      CMy = coder::interp2_dispatch(b_dv3, RudDefOut, psiT, rtNaN);
      CMz = coder::interp2_dispatch(b_dv4, RudDefOut, psiT, rtNaN);
    }
  }
  //  Transformation into UAV XYZ frame
  CFY_tmp = std::sin(vthr_z);
  b_CFY_tmp = std::cos(vthr_z);
  //  Aerodynamic Forces and Moments
  RudDefOut = ThrRPMOut / 60.0;
  RudDefOut = 1.225 * (RudDefOut * RudDefOut);
  MX_thr_tmp = RudDefOut * 0.001057227821024;
  MX_thr = MX_thr_tmp * CMx;
  //  Induced Velocity at Propeller Plane
  //  Induced velocity in hover
  d = 1.59 * (ThrRPMOut / 60.0) * 0.254;
  if ((vthr_x >= -0.2 * (0.5 * (d * 0.39538588745679831))) && (CFx > 0.0)) {
    //  if reverse velocity is greater than 20% of hover velocity
    psiT = 0.5 * (d * std::sqrt(CFx));
    MX_thr *= 0.2;
  } else {
    psiT = 0.0;
  }
  //  Thruster calculations done with Simulink blocks in original file
  //  Hacky Waqqas Khan patented division by two on Y and Z forces and moments
  RudDefOut *= 0.004162314256;
  F_Thr_idx_0 = RudDefOut * CFx;
  F_Thr_idx_1 = RudDefOut * (-CFy * b_CFY_tmp + 0.0 * CFY_tmp) * 0.5;
  F_Thr_idx_2 = RudDefOut * (-CFy * CFY_tmp - 0.0 * b_CFY_tmp) * 0.5;
  //  Gyroscopic moment
  W_Thr_SI = ThrRPMOut * 2.0 * 3.1415926535897931 / 60.0;
  //  Slipstream model
  //  Equations used from propeller slipstream paper
  //  Eq. 11
  //  Eq. 12
  vthr_z = 2.0 * psiT * 0.91823899371069173;
  //  Efflux velocity from Vi0_avg
  for (i = 0; i < 18; i++) {
    Cw_tmp = 7 * (i + 4);
    d = std::abs(Cts_Geom[Cw_tmp + 4]) * 0.001;
    d1 = Cts_Geom[Cw_tmp + 5];
    d2 = Cts_Geom[Cw_tmp + 6];
    d2 *= d2;
    VProp_Axial[i] = d2;
    d1 = std::sqrt(d1 * d1 + d2) * 0.001;
    if (d < 0.194056) {
      RudDefOut = d / 0.127;
      d2 = psiT * (d / 0.127 / std::sqrt(RudDefOut * RudDefOut + 1.0) + 1.0);
      VProp_Axial[i] = d2;
      //  Eq. 1
    } else if (d < 0.319532) {
      //  Eq. 13
      //  Eq. 14
      RudDefOut = (d1 - 0.058946599999999995 *
                            (1.0 - 0.1294 * (d - 0.194056) / 0.18796)) /
                  (0.1326 * ((d - 0.194056) - 0.09398) + 0.05210289974);
      d2 = vthr_z * (1.24 - 0.0765 * (d - 0.194056) / 0.18796) *
           std::exp(-(RudDefOut * RudDefOut));
      VProp_Axial[i] = d2;
      //  Eq. 15
    } else if (d < 0.79882999999999993) {
      //  Eq. 16
      //  Eq. 17
      RudDefOut = (d1 - 0.058946599999999995 *
                            (1.3 - 0.3059 * (d - 0.194056) / 0.18796)) /
                  (0.2295 * ((d - 0.194056) - 0.09398) + 0.030510760159999994);
      d2 = vthr_z * (1.37 - 0.1529 * (d - 0.194056) / 0.18796) *
           std::exp(-(RudDefOut * RudDefOut));
      VProp_Axial[i] = d2;
      //  Eq. 18
    } else {
      //  Eq. 19
      RudDefOut = d1 / (0.2411 * (d - 0.194056));
      d2 = vthr_z * (0.89 - 0.04 * (d - 0.194056) / 0.18796) *
           std::exp(-(RudDefOut * RudDefOut));
      VProp_Axial[i] = d2;
      //  Eq. 21
    }
    if (d2 < 0.01) {
      VProp_Axial[i] = 0.0;
    }
  }
  //  Calculate the radial component on Wing only
  //  DesR_Wing = DesR(1:6);
  //  VProp_Rad = zeros(size(DesR_Wing));
  //
  //  r = [0,0.35759,0.41593,...
  //      0.47616,0.5345,0.59849,0.65684,0.71518,0.77353,0.83375,0.89021]*Rp;
  //
  //  Vtip = Rp*RPM*2*pi/60;
  //  Vt = [0,0.035691875,...
  //      0.031344793,0.02951442,0.03203119,0.031344793,0.034319141,0.032946346,0.029285662,0.021964231,0.017845968]*Vtip;
  //
  //  for i = 1:length(DesR_Wing)
  //      if DesR_Wing(i) <= 0.89*Rp
  //          VProp_Rad(i) = interp1(r,Vt,DesR_Wing(i));
  //      else
  //          VProp_Rad(i) = 0;
  //      end
  //  end
  //  McFoamy Aerodynamics (Final forces and moments block)
  //  Density of air
  //  nu = 1.56e-05;  % Kinematic viscosity of air
  //  Re_crit = 3.e5;
  // --------------------------------------------------------------------------
  //                 Geometric properties of components
  // --------------------------------------------------------------------------
  //  AIRCRAFT CG
  //  Measured w.r.t the nose (from propeller plane)
  //  THRUSTER
  //  Thrust is assumed to act in x direction only
  //  WING GEOMETRY
  // Geometry(1,1);
  //  Aspect ratio
  //  TAIL GEOMETRY
  // Geometry(1,2);
  //  Aspect ratio
  Ct_idx_0 = Cts_Geom[79] * 0.001;
  Ctf_idx_0 = Cts_Geom[80] * 0.001;
  xt_idx_0 = Cts_Geom[81] * 0.001 - -0.27;
  d = Cts_Geom[82] * 0.001;
  yts_idx_0 = d;
  ytp_idx_0 = -d;
  zt_idx_0 = Cts_Geom[83] * 0.001 - 0.0058899999999999994;
  Ct_idx_1 = Cts_Geom[86] * 0.001;
  Ctf_idx_1 = Cts_Geom[87] * 0.001;
  xt_idx_1 = Cts_Geom[88] * 0.001 - -0.27;
  d = Cts_Geom[89] * 0.001;
  yts_idx_1 = d;
  ytp_idx_1 = -d;
  zt_idx_1 = Cts_Geom[90] * 0.001 - 0.0058899999999999994;
  Ct_idx_2 = Cts_Geom[93] * 0.001;
  Ctf_idx_2 = Cts_Geom[94] * 0.001;
  xt_idx_2 = Cts_Geom[95] * 0.001 - -0.27;
  d = Cts_Geom[96] * 0.001;
  zt_idx_2 = Cts_Geom[97] * 0.001 - 0.0058899999999999994;
  //  RUDDER GEOMETRY
  // Geometry(1,3);
  //  Aspect ratio
  //  BODY GEOMETRY
  // Geometry(1,4);
  //  Aspect ratio
  Cr_idx_0 = Cts_Geom[100] * 0.001;
  xr_idx_0 = Cts_Geom[102] * 0.001 - -0.27;
  yr_idx_0 = Cts_Geom[103] * 0.001;
  zr_idx_0 = Cts_Geom[104] * 0.001 - 0.0058899999999999994;
  CB_idx_0 = Cts_Geom[128] * 0.001;
  xB_idx_0 = Cts_Geom[130] * 0.001 - -0.27;
  yB_idx_0 = Cts_Geom[131] * 0.001;
  zB_idx_0 = Cts_Geom[132] * 0.001 - 0.0058899999999999994;
  Cr_idx_1 = Cts_Geom[107] * 0.001;
  xr_idx_1 = Cts_Geom[109] * 0.001 - -0.27;
  yr_idx_1 = Cts_Geom[110] * 0.001;
  zr_idx_1 = Cts_Geom[111] * 0.001 - 0.0058899999999999994;
  CB_idx_1 = Cts_Geom[135] * 0.001;
  xB_idx_1 = Cts_Geom[137] * 0.001 - -0.27;
  yB_idx_1 = Cts_Geom[138] * 0.001;
  zB_idx_1 = Cts_Geom[139] * 0.001 - 0.0058899999999999994;
  Cr_idx_2 = Cts_Geom[114] * 0.001;
  xr_idx_2 = Cts_Geom[116] * 0.001 - -0.27;
  yr_idx_2 = Cts_Geom[117] * 0.001;
  zr_idx_2 = Cts_Geom[118] * 0.001 - 0.0058899999999999994;
  CB_idx_2 = Cts_Geom[142] * 0.001;
  xB_idx_2 = Cts_Geom[144] * 0.001 - -0.27;
  yB_idx_2 = Cts_Geom[145] * 0.001;
  zB_idx_2 = Cts_Geom[146] * 0.001 - 0.0058899999999999994;
  Cr_idx_3 = Cts_Geom[121] * 0.001;
  xr_idx_3 = Cts_Geom[123] * 0.001 - -0.27;
  yr_idx_3 = Cts_Geom[124] * 0.001;
  zr_idx_3 = Cts_Geom[125] * 0.001 - 0.0058899999999999994;
  CB_idx_3 = Cts_Geom[149] * 0.001;
  xB_idx_3 = Cts_Geom[151] * 0.001 - -0.27;
  yB_idx_3 = Cts_Geom[152] * 0.001;
  zB_idx_3 = Cts_Geom[153] * 0.001 - 0.0058899999999999994;
  // --------------------------------------------------------------------------
  //                        Aerodynamic Constants
  // --------------------------------------------------------------------------
  //  alpStallP = 10*pi/180;
  //  alpStallN = -10*pi/180;
  //  HighAlpStart = 22*pi/180;
  //  HighAlpEnd = 158*pi/180;
  //  StallMdl = 'FullStall';
  // --------------------------------------------------------------------------
  //                        Wind disturbances
  // --------------------------------------------------------------------------
  //  Wind disturbance in body frame
  //  Vw_u = V_wind(1);
  //  Vw_v = V_wind(2);
  //  Vw_w = V_wind(3);
  //  % Vw_p = V_wind(4);
  //  % Vw_q = V_wind(5);
  //  % Vw_r = V_wind(6);
  //  Vw_p = 0;
  //  Vw_q = 0;
  //  Vw_r = 0;
  // --------------------------------------------------------------------------
  //      Body axis components of vehicle velocity relative to the air
  // --------------------------------------------------------------------------
  //  Assuming wind velocity is in body frame
  // --------------------------------------------------------------------------
  //               Slipstream velocities at reference points
  // --------------------------------------------------------------------------
  //  landing gears
  // --------------------------------------------------------------------------
  //                              Control Inputs
  // --------------------------------------------------------------------------
  RAilDef = 0.017453292519943295 * AilDefOut;
  AilDefOut = 0.017453292519943295 * ElevDefOut;
  // --------------------------------------------------------------------------
  //                    Wing, Tail, Rudder and Body velocities
  // --------------------------------------------------------------------------
  for (i = 0; i < 7; i++) {
    Cw_tmp = 7 * (i + 4);
    Cw[i] = Cts_Geom[Cw_tmp + 2] * 0.001;
    Cwf[i] = Cts_Geom[Cw_tmp + 3] * 0.001;
    d1 = Cts_Geom[Cw_tmp + 4] * 0.001 - -0.27;
    xw[i] = d1;
    d2 = Cts_Geom[Cw_tmp + 5] * 0.001;
    yws[i] = d2;
    ywp[i] = -d2;
    d3 = Cts_Geom[Cw_tmp + 6] * 0.001 - 0.0058899999999999994;
    zw[i] = d3;
    d3 = v_u + w_q * d3;
    d4 = VProp_Axial[i];
    d5 = (d3 - w_r * d2) + d4;
    vws_xcomp[i] = d5;
    d1 *= w_q;
    d6 = (v_w + w_p * d2) - d1;
    vws_zcomp[i] = d6;
    vws_xz[i] = std::sqrt(d5 * d5 + d6 * d6);
    d3 = (d3 - w_r * -d2) + d4;
    vwp_xcomp[i] = d3;
    d1 = (v_w + w_p * -d2) - d1;
    vwp_zcomp[i] = d1;
    vwp_xz[i] = std::sqrt(d3 * d3 + d1 * d1);
  }
  d1 = v_u + w_q * zt_idx_0;
  d2 = (d1 - w_r * yts_idx_0) + VProp_Axial[7];
  vts_xcomp_idx_0 = d2;
  d3 = w_q * xt_idx_0;
  d4 = (v_w + w_p * yts_idx_0) - d3;
  vts_zcomp_idx_0 = d4;
  ElevDefOut = std::sqrt(d2 * d2 + d4 * d4);
  d1 = (d1 - w_r * ytp_idx_0) + VProp_Axial[7];
  vtp_xcomp_idx_0 = d1;
  d3 = (v_w + w_p * ytp_idx_0) - d3;
  vtp_zcomp_idx_0 = d3;
  vtp_xz_idx_0 = std::sqrt(d1 * d1 + d3 * d3);
  d1 = v_u + w_q * zt_idx_1;
  d2 = (d1 - w_r * yts_idx_1) + VProp_Axial[8];
  vts_xcomp_idx_1 = d2;
  d3 = w_q * xt_idx_1;
  d4 = (v_w + w_p * yts_idx_1) - d3;
  vts_zcomp_idx_1 = d4;
  vts_xz_idx_1 = std::sqrt(d2 * d2 + d4 * d4);
  d1 = (d1 - w_r * ytp_idx_1) + VProp_Axial[8];
  vtp_xcomp_idx_1 = d1;
  d3 = (v_w + w_p * ytp_idx_1) - d3;
  vtp_zcomp_idx_1 = d3;
  vtp_xz_idx_1 = std::sqrt(d1 * d1 + d3 * d3);
  d1 = v_u + w_q * zt_idx_2;
  d2 = (d1 - w_r * d) + VProp_Axial[9];
  d3 = w_q * xt_idx_2;
  d4 = (v_w + w_p * d) - d3;
  vts_xz_idx_2 = std::sqrt(d2 * d2 + d4 * d4);
  d1 = (d1 - w_r * -d) + VProp_Axial[9];
  d3 = (v_w + w_p * -d) - d3;
  vtp_xz_idx_2 = std::sqrt(d1 * d1 + d3 * d3);
  RudDefOut = ((v_u + w_q * zr_idx_0) - w_r * yr_idx_0) + VProp_Axial[10];
  vr_ycomp[0] = zr_idx_0 + 0.0058899999999999994;
  psiT = ((v_u + w_q * zr_idx_1) - w_r * yr_idx_1) + VProp_Axial[11];
  vr_ycomp[1] = zr_idx_1 + 0.0058899999999999994;
  CFx = ((v_u + w_q * zr_idx_2) - w_r * yr_idx_2) + VProp_Axial[12];
  vr_ycomp[2] = zr_idx_2 + 0.0058899999999999994;
  CFy = ((v_u + w_q * zr_idx_3) - w_r * yr_idx_3) + VProp_Axial[13];
  vr_ycomp[3] = zr_idx_3 + 0.0058899999999999994;
  coder::b_sign(vr_ycomp);
  d5 = ((v_v + w_r * xr_idx_0) - w_p * zr_idx_0) + vr_ycomp[0] * 0.0;
  vr_ycomp[0] = d5;
  vr_xy_idx_0 = std::sqrt(RudDefOut * RudDefOut + d5 * d5);
  vthr_z = ((v_u + w_q * zB_idx_0) - w_r * yB_idx_0) + VProp_Axial[14];
  vB_ycomp[0] = zB_idx_0 + 0.0058899999999999994;
  d5 = ((v_v + w_r * xr_idx_1) - w_p * zr_idx_1) + vr_ycomp[1] * 0.0;
  vr_ycomp[1] = d5;
  vr_xy_idx_1 = std::sqrt(psiT * psiT + d5 * d5);
  CMx = ((v_u + w_q * zB_idx_1) - w_r * yB_idx_1) + VProp_Axial[15];
  vB_ycomp[1] = zB_idx_1 + 0.0058899999999999994;
  d5 = ((v_v + w_r * xr_idx_2) - w_p * zr_idx_2) + vr_ycomp[2] * 0.0;
  vr_ycomp[2] = d5;
  vr_xy_idx_2 = std::sqrt(CFx * CFx + d5 * d5);
  vthr_x = ((v_u + w_q * zB_idx_2) - w_r * yB_idx_2) + VProp_Axial[16];
  vB_ycomp[2] = zB_idx_2 + 0.0058899999999999994;
  d5 = ((v_v + w_r * xr_idx_3) - w_p * zr_idx_3) + vr_ycomp[3] * 0.0;
  vr_xy_idx_3 = std::sqrt(CFy * CFy + d5 * d5);
  ThrRPMOut = ((v_u + w_q * zB_idx_3) - w_r * yB_idx_3) + VProp_Axial[17];
  vB_ycomp[3] = zB_idx_3 + 0.0058899999999999994;
  coder::b_sign(vB_ycomp);
  d6 = ((v_v + w_r * xB_idx_0) - w_p * zB_idx_0) + vB_ycomp[0] * 0.0;
  vB_ycomp[0] = d6;
  vB_xy_idx_0 = std::sqrt(vthr_z * vthr_z + d6 * d6);
  d6 = ((v_v + w_r * xB_idx_1) - w_p * zB_idx_1) + vB_ycomp[1] * 0.0;
  vB_ycomp[1] = d6;
  vB_xy_idx_1 = std::sqrt(CMx * CMx + d6 * d6);
  d6 = ((v_v + w_r * xB_idx_2) - w_p * zB_idx_2) + vB_ycomp[2] * 0.0;
  vB_ycomp[2] = d6;
  vB_xy_idx_2 = std::sqrt(vthr_x * vthr_x + d6 * d6);
  d6 = ((v_v + w_r * xB_idx_3) - w_p * zB_idx_3) + vB_ycomp[3] * 0.0;
  vB_xy_idx_3 = std::sqrt(ThrRPMOut * ThrRPMOut + d6 * d6);
  // --------------------------------------------------------------------------
  //               Wing, Tail, Rudder and Body angles of attack
  // --------------------------------------------------------------------------
  //  Angle of attack = atan(V_zcomp / V_xcomp)
  //  range is -180 -> 180
  for (i = 0; i < 7; i++) {
    a_ws[i] = rt_atan2d_snf(vws_zcomp[i], vws_xcomp[i]);
    a_wp[i] = rt_atan2d_snf(vwp_zcomp[i], vwp_xcomp[i]);
  }
  double Crf_idx_2;
  double a_B_idx_1;
  double a_B_idx_2;
  double a_B_idx_3;
  double a_r_idx_1;
  double a_r_idx_2;
  double a_r_idx_3;
  a_ts_idx_0 = rt_atan2d_snf(vts_zcomp_idx_0, vts_xcomp_idx_0);
  a_tp_idx_0 = rt_atan2d_snf(vtp_zcomp_idx_0, vtp_xcomp_idx_0);
  a_ts_idx_1 = rt_atan2d_snf(vts_zcomp_idx_1, vts_xcomp_idx_1);
  a_tp_idx_1 = rt_atan2d_snf(vtp_zcomp_idx_1, vtp_xcomp_idx_1);
  a_ts_idx_2 = rt_atan2d_snf(d4, d2);
  a_tp_idx_2 = rt_atan2d_snf(d3, d1);
  //  Vertical AoA for rudder and body = atan(V_ycomp / V_xcomp)
  //  range is -180 -> 180
  a_r_idx_0 = rt_atan2d_snf(vr_ycomp[0], RudDefOut);
  a_B_idx_0 = rt_atan2d_snf(vB_ycomp[0], vthr_z);
  a_r_idx_1 = rt_atan2d_snf(vr_ycomp[1], psiT);
  a_B_idx_1 = rt_atan2d_snf(vB_ycomp[1], CMx);
  a_r_idx_2 = rt_atan2d_snf(vr_ycomp[2], CFx);
  a_B_idx_2 = rt_atan2d_snf(vB_ycomp[2], vthr_x);
  a_r_idx_3 = rt_atan2d_snf(d5, CFy);
  a_B_idx_3 = rt_atan2d_snf(d6, ThrRPMOut);
  // --------------------------------------------------------------------------
  //                    Wing lift and drag coefficient
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------
  //                        Wing-Tail interference
  // --------------------------------------------------------------------------
  //  % Bound vortices of wing
  //  G_ws = 0.5*Cw.*CL_ws.*vws_xz;
  //  G_wp = 0.5*Cw.*CL_wp.*vwp_xz;
  //
  //  % Induced downwash on tail due to bound vortices
  //  ViTS_ws = zeros(Nt,Nw);
  //  ViTS_wp = zeros(Nt,Nw);
  //  ViTP_ws = zeros(Nt,Nw);
  //  ViTP_wp = zeros(Nt,Nw);
  //
  //  for i = 1:Nt
  //      for j = 1:Nw
  //          ViTS_ws(i,j) = G_ws(j)*GTermTS_ws(i,j);
  //          ViTS_wp(i,j) = G_wp(j)*GTermTS_wp(i,j);
  //          ViTP_ws(i,j) = G_ws(j)*GTermTP_ws(i,j);
  //          ViTP_wp(i,j) = G_wp(j)*GTermTP_wp(i,j);
  //      end
  //  end
  //
  //  % Trailing vortices of wing
  //  delG_ws = zeros(1,Nw+2);
  //  delG_wp = zeros(1,Nw+2);
  //
  //  delG_ws(1) = G_ws(1) - 0;
  //  delG_wp(1) = 0 - G_wp(1);
  //
  //  for i = 2:Nw-1
  //      delG_ws(i) = G_ws(i) - G_ws(i-1);
  //      delG_wp(i) = G_wp(i-1) - G_wp(i);
  //  end
  //
  //  delG_ws(Nw) = G_ws(Nw);
  //  delG_ws(Nw+1) = 0 - G_ws(Nw-1);
  //  delG_ws(Nw+2) = 0 - G_ws(Nw);
  //
  //  delG_wp(Nw) = 0 - G_wp(Nw);
  //  delG_wp(Nw+1) = G_wp(Nw-1) - 0;
  //  delG_wp(Nw+2) = G_wp(Nw) - 0;
  //
  //  % Induced downwash on tail due to trailing vortices
  //  delViTS_ws = zeros(Nt,Nw+2);
  //  delViTS_wp = zeros(Nt,Nw+2);
  //  delViTP_ws = zeros(Nt,Nw+2);
  //  delViTP_wp = zeros(Nt,Nw+2);
  //
  //  for i = 1:Nt
  //      for j = 1:Nw+2
  //          delViTS_ws(i,j) = delG_ws(j)*delGTermTS_ws(i,j);
  //          delViTS_wp(i,j) = delG_wp(j)*delGTermTS_wp(i,j);
  //          delViTP_ws(i,j) = delG_ws(j)*delGTermTP_ws(i,j);
  //          delViTP_wp(i,j) = delG_wp(j)*delGTermTP_wp(i,j);
  //      end
  //  end
  //
  //  % Summing all induced downwash
  //  ViTS = zeros(1,Nt);
  //  ViTP = zeros(1,Nt);
  //
  //  for i = 1:Nt
  //      ViTS(i) = sum(ViTS_ws(i,:)) + sum(ViTS_wp(i,:)) + sum(delViTS_ws(i,:))
  //      + sum(delViTS_wp(i,:)); ViTP(i) = sum(ViTP_ws(i,:)) +
  //      sum(ViTP_wp(i,:)) + sum(delViTP_ws(i,:)) + sum(delViTP_wp(i,:));
  //  end
  //
  //  % Induced angle of attack on tail
  //  vts_zcomp = vts_zcomp - ViTS;
  //  vtp_zcomp = vtp_zcomp - ViTP;
  //
  //  a_ts = atan2(vts_zcomp,vts_xcomp);
  //  a_tp = atan2(vtp_zcomp,vtp_xcomp);
  // --------------------------------------------------------------------------
  //                   Tail's lift and drag coefficient
  // --------------------------------------------------------------------------
  b_flappedAirfoil_LAR_v1_FlatPlate(a_ts_idx_0, Ctf_idx_0, Ct_idx_0, AilDefOut,
                                    1.98, &d1, &d2, &d3, &d4);
  CL_ts_idx_0 = d2;
  vtp_xcomp_idx_0 = d3;
  CM_ts_idx_0 = d4;
  b_flappedAirfoil_LAR_v1_FlatPlate(a_tp_idx_0, Ctf_idx_0, Ct_idx_0, AilDefOut,
                                    1.98, &d1, &d5, &d6, &d7);
  CL_tp_idx_0 = d5;
  vtp_zcomp_idx_0 = d6;
  CM_tp_idx_0 = d7;
  b_flappedAirfoil_LAR_v1_FlatPlate(a_ts_idx_1, Ctf_idx_1, Ct_idx_1, AilDefOut,
                                    1.98, &d1, &d2, &d3, &d4);
  CL_ts_idx_1 = d2;
  vtp_xcomp_idx_1 = d3;
  CM_ts_idx_1 = d4;
  b_flappedAirfoil_LAR_v1_FlatPlate(a_tp_idx_1, Ctf_idx_1, Ct_idx_1, AilDefOut,
                                    1.98, &d1, &d5, &d6, &d7);
  CL_tp_idx_1 = d5;
  vtp_zcomp_idx_1 = d6;
  CM_tp_idx_1 = d7;
  b_flappedAirfoil_LAR_v1_FlatPlate(a_ts_idx_2, Ctf_idx_2, Ct_idx_2, AilDefOut,
                                    1.98, &d1, &d2, &d3, &d4);
  b_flappedAirfoil_LAR_v1_FlatPlate(a_tp_idx_2, Ctf_idx_2, Ct_idx_2, AilDefOut,
                                    1.98, &d1, &d5, &d6, &d7);
  // --------------------------------------------------------------------------
  //                        Wing-Rudder interference
  // --------------------------------------------------------------------------
  //  % Induced downwash on rudder due to bound vortices
  //  ViR_ws_xcomp = zeros(Nr,Nw);
  //  ViR_ws_zcomp = zeros(Nr,Nw);
  //  ViR_wp_xcomp = zeros(Nr,Nw);
  //  ViR_wp_zcomp = zeros(Nr,Nw);
  //
  //  for i = 1:Nr
  //      for j = 1:Nw
  //          ViR_ws_xcomp(i,j) = G_ws(j)*GTermR_w_xcomp(i,j);
  //          ViR_ws_zcomp(i,j) = G_ws(j)*GTermR_w_zcomp(i,j);
  //
  //          ViR_wp_xcomp(i,j) = G_wp(j)*GTermR_w_xcomp(i,j);
  //          ViR_wp_zcomp(i,j) = G_wp(j)*GTermR_w_zcomp(i,j);
  //      end
  //  end
  //
  //  % Induced downwash on rudder due to trailing vortices
  //  delViR_ws_ycomp = zeros(Nr,Nw+1);
  //  delViR_ws_zcomp = zeros(Nr,Nw+1);
  //  delViR_wp_ycomp = zeros(Nr,Nw+1);
  //  delViR_wp_zcomp = zeros(Nr,Nw+1);
  //
  //  for i = 1:Nr
  //      for j = 1:Nw+1
  //          delViR_ws_ycomp(i,j) = delG_ws(j)*delGTermR_ws_ycomp(i,j);
  //          delViR_ws_zcomp(i,j) = delG_ws(j)*delGTermR_ws_zcomp(i,j);
  //
  //          delViR_wp_ycomp(i,j) = delG_wp(j)*delGTermR_wp_ycomp(i,j);
  //          delViR_wp_zcomp(i,j) = delG_wp(j)*delGTermR_wp_zcomp(i,j);
  //      end
  //  end
  //
  //  % Summing all induced downwash
  //  ViR_xcomp = zeros(1,Nr);
  //  ViR_ycomp = zeros(1,Nr);
  //  ViR_zcomp = zeros(1,Nr);
  //
  //  for i = 1:Nr
  //      ViR_xcomp(i) = sum(ViR_ws_xcomp(i,:)) + sum(ViR_wp_xcomp(i,:));
  //      ViR_ycomp(i) = sum(delViR_ws_ycomp(i,:)) + sum(delViR_wp_ycomp(i,:));
  //      ViR_zcomp(i) = sum(ViR_ws_zcomp(i,:)) + sum(ViR_wp_zcomp(i,:)) +
  //      sum(delViR_ws_zcomp(i,:)) + sum(delViR_wp_zcomp(i,:));
  //  end
  //
  //  % Induced angle of attack on rudder
  //  vr_xcomp = vr_xcomp - ViR_xcomp;
  //  vr_ycomp = vr_ycomp - ViR_ycomp;
  //  vr_zcomp = vr_zcomp - ViR_zcomp;
  //
  //  a_r = atan2(vr_ycomp,vr_xcomp);
  // --------------------------------------------------------------------------
  //                   Rudder lift and drag coefficient
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------
  //                   Fuselage lift and drag coefficient
  // --------------------------------------------------------------------------
  //  check sign of deflection into flapped code
  c_flappedAirfoil_LAR_v1_FlatPlate(a_r_idx_0, Cts_Geom[101] * 0.001, Cr_idx_0,
                                    -CtrlDef_idx_3, 1.98, &d1, &d8, &d9, &d10);
  CL_r_idx_0 = d8;
  RudDefOut = d9;
  CM_r_idx_0 = d10;
  flappedAirfoil_LAR_v1_FlatPlate(a_B_idx_0, Cts_Geom[129] * 0.001, CB_idx_0,
                                  1.98, &d1, &d11, &d12, &d13);
  CL_B_idx_0 = d11;
  vthr_z = d12;
  CM_B_idx_0 = d13;
  //  check sign of deflection into flapped code
  c_flappedAirfoil_LAR_v1_FlatPlate(a_r_idx_1, Cts_Geom[108] * 0.001, Cr_idx_1,
                                    -CtrlDef_idx_3, 1.98, &d1, &d8, &d9, &d10);
  CL_r_idx_1 = d8;
  psiT = d9;
  CM_r_idx_1 = d10;
  flappedAirfoil_LAR_v1_FlatPlate(a_B_idx_1, Cts_Geom[136] * 0.001, CB_idx_1,
                                  1.98, &d1, &d11, &d12, &d13);
  CL_B_idx_1 = d11;
  ThrRPMOut = d12;
  CM_B_idx_1 = d13;
  //  check sign of deflection into flapped code
  c_flappedAirfoil_LAR_v1_FlatPlate(a_r_idx_2, Cts_Geom[115] * 0.001, Cr_idx_2,
                                    -CtrlDef_idx_3, 1.98, &d1, &d8, &d9, &d10);
  CL_r_idx_2 = d8;
  AilDefOut = d9;
  CM_r_idx_2 = d10;
  flappedAirfoil_LAR_v1_FlatPlate(a_B_idx_2, Cts_Geom[143] * 0.001, CB_idx_2,
                                  1.98, &d1, &d11, &d12, &d13);
  CL_B_idx_2 = d11;
  Crf_idx_2 = d12;
  CM_B_idx_2 = d13;
  //  check sign of deflection into flapped code
  c_flappedAirfoil_LAR_v1_FlatPlate(a_r_idx_3, Cts_Geom[122] * 0.001, Cr_idx_3,
                                    -CtrlDef_idx_3, 1.98, &d1, &d8, &d9, &d10);
  flappedAirfoil_LAR_v1_FlatPlate(a_B_idx_3, Cts_Geom[150] * 0.001, CB_idx_3,
                                  1.98, &d1, &d11, &d12, &d13);
  // --------------------------------------------------------------------------
  //                   Structural rods
  // --------------------------------------------------------------------------
  //  tol=1e-6;
  //  omega=[p,q,r]';
  //  V=[u,v,w]';
  //  CD_lg=1.1;
  //  D_rods=zeros(3,1);
  //  M_rods=zeros(3,1);
  //  for i=1:size(rods,1)
  //      ai=rods(i,1:3)';  % pt1
  //      bi=rods(i,4:6)';  % pt2
  //      di=rods(i,7);     % diameter of the rod
  //      ri=(ai+bi)/2;   % center of the rod w/r to the CG
  //      li=bi-ai;       % length vector of the rod
  //      Vi=V+cross(omega,ri)+[Vp_lg(i),0,0]';
  //      sin_th=norm(cross(Vi,li))/max(norm(Vi),tol)/norm(li);
  //      V_l_l=cross(cross(Vi,li),li);
  //      if Vi'*V_l_l<0  % Vi pointing x forward, Nui pointing backward
  //          Nui=V_l_l/max(norm(Vi)*norm(li)^2,tol);
  //      else
  //          Nui=-V_l_l/max(norm(Vi)*norm(li)^2,tol);
  //      end
  //      Di=0.5*rho*norm(Vi)^2*norm(li)*di*CD_lg*sin_th^2*Nui;
  //      Mi=cross(ri,Di);
  //      D_rods=D_rods+Di;
  //      M_rods=M_rods+Mi;
  //  end
  // --------------------------------------------------------------------------
  //                    Aerodynamic Forces & Moments
  // --------------------------------------------------------------------------
  //  x-direction forces
  Fx_ts_tmp_idx_2 = 0.6125 * (Cts_Geom[77] * 0.001) * Ct_idx_0;
  d1 = std::cos(a_ts_idx_0);
  vts_zcomp_idx_0 = d1;
  d14 = std::sin(a_ts_idx_0);
  a_ts_idx_0 = d14;
  y_idx_2 = Fx_ts_tmp_idx_2 * (ElevDefOut * ElevDefOut);
  y_idx_0 = y_idx_2;
  Ctf_idx_0 = y_idx_2 * (CL_ts_idx_0 * d14 - vtp_xcomp_idx_0 * d1);
  d15 = std::cos(a_tp_idx_0);
  CtrlDef_idx_3 = d15;
  d16 = std::sin(a_tp_idx_0);
  a_tp_idx_0 = d16;
  Fx_ts_tmp_idx_2 *= vtp_xz_idx_0 * vtp_xz_idx_0;
  Fx_ts_tmp_idx_0 = Fx_ts_tmp_idx_2;
  vts_xcomp_idx_0 =
      Fx_ts_tmp_idx_2 * (CL_tp_idx_0 * d16 - vtp_zcomp_idx_0 * d15);
  Fx_ts_tmp_idx_2 = 0.6125 * (Cts_Geom[84] * 0.001) * Ct_idx_1;
  d1 = std::cos(a_ts_idx_1);
  vts_zcomp_idx_1 = d1;
  d14 = std::sin(a_ts_idx_1);
  a_ts_idx_1 = d14;
  y_idx_2 = Fx_ts_tmp_idx_2 * (vts_xz_idx_1 * vts_xz_idx_1);
  y_idx_1 = y_idx_2;
  Ctf_idx_1 = y_idx_2 * (CL_ts_idx_1 * d14 - vtp_xcomp_idx_1 * d1);
  d15 = std::cos(a_tp_idx_1);
  ElevDefOut = d15;
  d16 = std::sin(a_tp_idx_1);
  a_tp_idx_1 = d16;
  Fx_ts_tmp_idx_2 *= vtp_xz_idx_1 * vtp_xz_idx_1;
  Fx_ts_tmp_idx_1 = Fx_ts_tmp_idx_2;
  vts_xcomp_idx_1 =
      Fx_ts_tmp_idx_2 * (CL_tp_idx_1 * d16 - vtp_zcomp_idx_1 * d15);
  Fx_ts_tmp_idx_2 = 0.6125 * (Cts_Geom[91] * 0.001) * Ct_idx_2;
  d1 = std::cos(a_ts_idx_2);
  d14 = std::sin(a_ts_idx_2);
  y_idx_2 = Fx_ts_tmp_idx_2 * (vts_xz_idx_2 * vts_xz_idx_2);
  Ctf_idx_2 = y_idx_2 * (d2 * d14 - d3 * d1);
  d15 = std::cos(a_tp_idx_2);
  d16 = std::sin(a_tp_idx_2);
  Fx_ts_tmp_idx_2 *= vtp_xz_idx_2 * vtp_xz_idx_2;
  vts_xcomp_idx_2 = Fx_ts_tmp_idx_2 * (d5 * d16 - d6 * d15);
  //  y-direction forces
  //  INCLUDE FRICTION DRAG
  //  INCLUDE FRICTION DRAG
  //  INCLUDE FRICTION DRAG
  //  INCLUDE FRICTION DRAG
  CFy = std::cos(a_r_idx_0);
  vthr_x = std::sin(a_r_idx_0);
  y_idx_3 =
      0.6125 * (Cts_Geom[98] * 0.001) * Cr_idx_0 * (vr_xy_idx_0 * vr_xy_idx_0);
  b_y_idx_0 = y_idx_3;
  vB_ycomp[0] = y_idx_3 * (CL_r_idx_0 * vthr_x - RudDefOut * CFy);
  CFx = std::cos(a_B_idx_0);
  CMx = std::sin(a_B_idx_0);
  vr_xy_idx_0 =
      0.6125 * (Cts_Geom[126] * 0.001) * CB_idx_0 * (vB_xy_idx_0 * vB_xy_idx_0);
  a_B_idx_0 = vr_xy_idx_0;
  vr_ycomp[0] = vr_xy_idx_0 * (CL_B_idx_0 * CMx - vthr_z * CFx);
  CL_r_idx_0 = y_idx_3 * (-CL_r_idx_0 * CFy - RudDefOut * vthr_x);
  CL_B_idx_0 = vr_xy_idx_0 * (-CL_B_idx_0 * CFx - vthr_z * CMx);
  CFy = std::cos(a_r_idx_1);
  vthr_x = std::sin(a_r_idx_1);
  y_idx_3 =
      0.6125 * (Cts_Geom[105] * 0.001) * Cr_idx_1 * (vr_xy_idx_1 * vr_xy_idx_1);
  a_r_idx_0 = y_idx_3;
  vB_ycomp[1] = y_idx_3 * (CL_r_idx_1 * vthr_x - psiT * CFy);
  CFx = std::cos(a_B_idx_1);
  CMx = std::sin(a_B_idx_1);
  vr_xy_idx_0 =
      0.6125 * (Cts_Geom[133] * 0.001) * CB_idx_1 * (vB_xy_idx_1 * vB_xy_idx_1);
  a_tp_idx_2 = vr_xy_idx_0;
  vr_ycomp[1] = vr_xy_idx_0 * (CL_B_idx_1 * CMx - ThrRPMOut * CFx);
  CL_r_idx_1 = y_idx_3 * (-CL_r_idx_1 * CFy - psiT * vthr_x);
  CL_B_idx_1 = vr_xy_idx_0 * (-CL_B_idx_1 * CFx - ThrRPMOut * CMx);
  CFy = std::cos(a_r_idx_2);
  vthr_x = std::sin(a_r_idx_2);
  y_idx_3 =
      0.6125 * (Cts_Geom[112] * 0.001) * Cr_idx_2 * (vr_xy_idx_2 * vr_xy_idx_2);
  a_ts_idx_2 = y_idx_3;
  vB_ycomp[2] = y_idx_3 * (CL_r_idx_2 * vthr_x - AilDefOut * CFy);
  CFx = std::cos(a_B_idx_2);
  CMx = std::sin(a_B_idx_2);
  vr_xy_idx_0 =
      0.6125 * (Cts_Geom[140] * 0.001) * CB_idx_2 * (vB_xy_idx_2 * vB_xy_idx_2);
  vtp_xz_idx_2 = vr_xy_idx_0;
  vr_ycomp[2] = vr_xy_idx_0 * (CL_B_idx_2 * CMx - Crf_idx_2 * CFx);
  CL_r_idx_2 = y_idx_3 * (-CL_r_idx_2 * CFy - AilDefOut * vthr_x);
  CL_B_idx_2 = vr_xy_idx_0 * (-CL_B_idx_2 * CFx - Crf_idx_2 * CMx);
  CFy = std::cos(a_r_idx_3);
  vthr_x = std::sin(a_r_idx_3);
  y_idx_3 =
      0.6125 * (Cts_Geom[119] * 0.001) * Cr_idx_3 * (vr_xy_idx_3 * vr_xy_idx_3);
  vB_ycomp[3] = y_idx_3 * (d8 * vthr_x - d9 * CFy);
  CFx = std::cos(a_B_idx_3);
  CMx = std::sin(a_B_idx_3);
  vr_xy_idx_0 =
      0.6125 * (Cts_Geom[147] * 0.001) * CB_idx_3 * (vB_xy_idx_3 * vB_xy_idx_3);
  vr_ycomp[3] = vr_xy_idx_0 * (d11 * CMx - d12 * CFx);
  vtp_xz_idx_1 = y_idx_3 * (-d8 * CFy - d9 * vthr_x);
  vts_xz_idx_2 = vr_xy_idx_0 * (-d11 * CFx - d12 * CMx);
  //  z-direction forces
  for (i = 0; i < 7; i++) {
    d8 = a_ws[i];
    d9 = Cwf[i];
    d11 = Cw[i];
    flappedAirfoil_LAR_v1_FlatPlate(d8, d9, d11, RAilDef, 1.98, &d12, &CFy,
                                    &vthr_x, &CFx);
    CM_ws[i] = CFx;
    d12 = a_wp[i];
    flappedAirfoil_LAR_v1_FlatPlate(d12, d9, d11, -RAilDef, 1.98, &CFx, &CMx,
                                    &RudDefOut, &vthr_z);
    CM_wp[i] = vthr_z;
    d9 = 0.6125 * (Cts_Geom[7 * (i + 4)] * 0.001) * d11;
    d11 = std::cos(d8);
    d8 = std::sin(d8);
    a_ws[i] = d8;
    CFx = vws_xz[i];
    CFx = d9 * (CFx * CFx);
    vwp_xcomp[i] = CFx;
    vwp_zcomp[i] = CFx * (CFy * d8 - vthr_x * d11);
    vthr_z = std::cos(d12);
    d12 = std::sin(d12);
    a_wp[i] = d12;
    psiT = vwp_xz[i];
    d9 *= psiT * psiT;
    Cwf[i] = d9;
    vws_xcomp[i] = d9 * (CMx * d12 - RudDefOut * vthr_z);
    CL_ws[i] = CFx * (-CFy * d11 - vthr_x * d8);
    CL_wp[i] = d9 * (-CMx * vthr_z - RudDefOut * d12);
  }
  CL_ts_idx_0 =
      y_idx_0 * (-CL_ts_idx_0 * vts_zcomp_idx_0 - vtp_xcomp_idx_0 * a_ts_idx_0);
  CL_tp_idx_0 = Fx_ts_tmp_idx_0 *
                (-CL_tp_idx_0 * CtrlDef_idx_3 - vtp_zcomp_idx_0 * a_tp_idx_0);
  CL_ts_idx_1 =
      y_idx_1 * (-CL_ts_idx_1 * vts_zcomp_idx_1 - vtp_xcomp_idx_1 * a_ts_idx_1);
  CL_tp_idx_1 = Fx_ts_tmp_idx_1 *
                (-CL_tp_idx_1 * ElevDefOut - vtp_zcomp_idx_1 * a_tp_idx_1);
  vts_xz_idx_1 = y_idx_2 * (-d2 * d1 - d3 * d14);
  RudDefOut = Fx_ts_tmp_idx_2 * (-d5 * d15 - d6 * d16);
  //  INCLUDE FRICTION DRAG
  //  INCLUDE FRICTION DRAG
  //  x-direction moments
  //  y-direction moments
  //  z-direction moments
  //  WHY NEGATIVE SIGN???
  //  WHY NEGATIVE SIGN???
  //  Total forces and moments in x, y and z directions
  //  Fx = F_Thr(1);
  //  Fy = F_Thr(2);
  //  Fz = F_Thr(3);
  //  F = [Fx,Fy,Fz]+D_rods';
  //  M = [Mx, My, Mz]+M_rods';
  CMx = vwp_zcomp[0];
  CFx = vws_xcomp[0];
  vthr_x = CL_ws[0];
  CFy = CL_wp[0];
  for (i = 0; i < 6; i++) {
    CMx += vwp_zcomp[i + 1];
    CFx += vws_xcomp[i + 1];
    vthr_x += CL_ws[i + 1];
    CFy += CL_wp[i + 1];
  }
  *Fx =
      (((((CMx + CFx) + ((Ctf_idx_0 + Ctf_idx_1) + Ctf_idx_2)) +
         ((vts_xcomp_idx_0 + vts_xcomp_idx_1) + vts_xcomp_idx_2)) +
        (((vB_ycomp[0] + vB_ycomp[1]) + vB_ycomp[2]) + vB_ycomp[3])) +
       (((vr_ycomp[0] + vr_ycomp[1]) + vr_ycomp[2]) + vr_ycomp[3])) +
      F_Thr_idx_0;
  *Fy =
      ((((CL_r_idx_0 + CL_r_idx_1) + CL_r_idx_2) + vtp_xz_idx_1) +
       (((CL_B_idx_0 + CL_B_idx_1) + CL_B_idx_2) + vts_xz_idx_2)) +
      F_Thr_idx_1;
  *Fz =
      (((vthr_x + CFy) + ((CL_ts_idx_0 + CL_ts_idx_1) + vts_xz_idx_1)) +
       ((CL_tp_idx_0 + CL_tp_idx_1) + RudDefOut)) +
      F_Thr_idx_2;
  for (Cw_tmp = 0; Cw_tmp < 7; Cw_tmp++) {
    vws_zcomp[Cw_tmp] = yws[Cw_tmp] * CL_ws[Cw_tmp] - zw[Cw_tmp] * 0.0;
  }
  CMx = vws_zcomp[0];
  for (i = 0; i < 6; i++) {
    CMx += vws_zcomp[i + 1];
  }
  for (Cw_tmp = 0; Cw_tmp < 7; Cw_tmp++) {
    vws_zcomp[Cw_tmp] = ywp[Cw_tmp] * CL_wp[Cw_tmp] - zw[Cw_tmp] * 0.0;
  }
  CFx = vws_zcomp[0];
  for (i = 0; i < 6; i++) {
    CFx += vws_zcomp[i + 1];
  }
  vthr_x = ((yts_idx_0 * CL_ts_idx_0 - zt_idx_0 * 0.0) +
            (yts_idx_1 * CL_ts_idx_1 - zt_idx_1 * 0.0)) +
           (d * vts_xz_idx_1 - zt_idx_2 * 0.0);
  CFy = ((ytp_idx_0 * CL_tp_idx_0 - zt_idx_0 * 0.0) +
         (ytp_idx_1 * CL_tp_idx_1 - zt_idx_1 * 0.0)) +
        (-d * RudDefOut - zt_idx_2 * 0.0);
  ThrRPMOut = (((yr_idx_0 * 0.0 - zr_idx_0 * CL_r_idx_0) +
                (yr_idx_1 * 0.0 - zr_idx_1 * CL_r_idx_1)) +
               (yr_idx_2 * 0.0 - zr_idx_2 * CL_r_idx_2)) +
              (yr_idx_3 * 0.0 - zr_idx_3 * vtp_xz_idx_1);
  AilDefOut = (((yB_idx_0 * 0.0 - zB_idx_0 * CL_B_idx_0) +
                (yB_idx_1 * 0.0 - zB_idx_1 * CL_B_idx_1)) +
               (yB_idx_2 * 0.0 - zB_idx_2 * CL_B_idx_2)) +
              (yB_idx_3 * 0.0 - zB_idx_3 * vts_xz_idx_2);
  for (Cw_tmp = 0; Cw_tmp < 7; Cw_tmp++) {
    CL_ws[Cw_tmp] =
        (zw[Cw_tmp] * vwp_zcomp[Cw_tmp] - xw[Cw_tmp] * CL_ws[Cw_tmp]) +
        vwp_xcomp[Cw_tmp] * Cw[Cw_tmp] * CM_ws[Cw_tmp];
  }
  ElevDefOut = CL_ws[0];
  for (i = 0; i < 6; i++) {
    ElevDefOut += CL_ws[i + 1];
  }
  for (Cw_tmp = 0; Cw_tmp < 7; Cw_tmp++) {
    zw[Cw_tmp] = (zw[Cw_tmp] * vws_xcomp[Cw_tmp] - xw[Cw_tmp] * CL_wp[Cw_tmp]) +
                 Cwf[Cw_tmp] * Cw[Cw_tmp] * CM_wp[Cw_tmp];
  }
  vtp_xz_idx_0 = zw[0];
  for (i = 0; i < 6; i++) {
    vtp_xz_idx_0 += zw[i + 1];
  }
  CL_ts_idx_0 = (zt_idx_0 * Ctf_idx_0 - xt_idx_0 * CL_ts_idx_0) +
                y_idx_0 * Ct_idx_0 * CM_ts_idx_0;
  CL_ts_idx_1 = (zt_idx_1 * Ctf_idx_1 - xt_idx_1 * CL_ts_idx_1) +
                y_idx_1 * Ct_idx_1 * CM_ts_idx_1;
  vts_xz_idx_1 = (zt_idx_2 * Ctf_idx_2 - xt_idx_2 * vts_xz_idx_1) +
                 y_idx_2 * Ct_idx_2 * d4;
  zt_idx_0 = (zt_idx_0 * vts_xcomp_idx_0 - xt_idx_0 * CL_tp_idx_0) +
             Fx_ts_tmp_idx_0 * Ct_idx_0 * CM_tp_idx_0;
  zt_idx_1 = (zt_idx_1 * vts_xcomp_idx_1 - xt_idx_1 * CL_tp_idx_1) +
             Fx_ts_tmp_idx_1 * Ct_idx_1 * CM_tp_idx_1;
  zt_idx_2 = (zt_idx_2 * vts_xcomp_idx_2 - xt_idx_2 * RudDefOut) +
             Fx_ts_tmp_idx_2 * Ct_idx_2 * d7;
  zr_idx_0 = zr_idx_0 * vB_ycomp[0] - xr_idx_0 * 0.0;
  zr_idx_1 = zr_idx_1 * vB_ycomp[1] - xr_idx_1 * 0.0;
  zr_idx_2 = zr_idx_2 * vB_ycomp[2] - xr_idx_2 * 0.0;
  zr_idx_3 = zr_idx_3 * vB_ycomp[3] - xr_idx_3 * 0.0;
  zB_idx_0 = zB_idx_0 * vr_ycomp[0] - xB_idx_0 * 0.0;
  zB_idx_1 = zB_idx_1 * vr_ycomp[1] - xB_idx_1 * 0.0;
  zB_idx_2 = zB_idx_2 * vr_ycomp[2] - xB_idx_2 * 0.0;
  zB_idx_3 = zB_idx_3 * vr_ycomp[3] - xB_idx_3 * 0.0;
  for (Cw_tmp = 0; Cw_tmp < 7; Cw_tmp++) {
    yws[Cw_tmp] = xw[Cw_tmp] * 0.0 - yws[Cw_tmp] * vwp_zcomp[Cw_tmp];
  }
  vthr_z = yws[0];
  for (i = 0; i < 6; i++) {
    vthr_z += yws[i + 1];
  }
  for (Cw_tmp = 0; Cw_tmp < 7; Cw_tmp++) {
    xw[Cw_tmp] = xw[Cw_tmp] * 0.0 - ywp[Cw_tmp] * vws_xcomp[Cw_tmp];
  }
  psiT = xw[0];
  for (i = 0; i < 6; i++) {
    psiT += xw[i + 1];
  }
  yts_idx_0 = xt_idx_0 * 0.0 - yts_idx_0 * Ctf_idx_0;
  yts_idx_1 = xt_idx_1 * 0.0 - yts_idx_1 * Ctf_idx_1;
  RudDefOut = (yts_idx_0 + yts_idx_1) + (xt_idx_2 * 0.0 - d * Ctf_idx_2);
  xt_idx_0 = xt_idx_0 * 0.0 - ytp_idx_0 * vts_xcomp_idx_0;
  xt_idx_1 = xt_idx_1 * 0.0 - ytp_idx_1 * vts_xcomp_idx_1;
  xt_idx_2 = xt_idx_2 * 0.0 - -d * vts_xcomp_idx_2;
  xr_idx_0 = (xr_idx_0 * CL_r_idx_0 - yr_idx_0 * vB_ycomp[0]) -
             b_y_idx_0 * Cr_idx_0 * CM_r_idx_0;
  xr_idx_1 = (xr_idx_1 * CL_r_idx_1 - yr_idx_1 * vB_ycomp[1]) -
             a_r_idx_0 * Cr_idx_1 * CM_r_idx_1;
  xr_idx_2 = (xr_idx_2 * CL_r_idx_2 - yr_idx_2 * vB_ycomp[2]) -
             a_ts_idx_2 * Cr_idx_2 * CM_r_idx_2;
  xr_idx_3 = (xr_idx_3 * vtp_xz_idx_1 - yr_idx_3 * vB_ycomp[3]) -
             y_idx_3 * Cr_idx_3 * d10;
  xB_idx_0 = (xB_idx_0 * CL_B_idx_0 - yB_idx_0 * vr_ycomp[0]) -
             a_B_idx_0 * CB_idx_0 * CM_B_idx_0;
  xB_idx_1 = (xB_idx_1 * CL_B_idx_1 - yB_idx_1 * vr_ycomp[1]) -
             a_tp_idx_2 * CB_idx_1 * CM_B_idx_1;
  xB_idx_2 = (xB_idx_2 * CL_B_idx_2 - yB_idx_2 * vr_ycomp[2]) -
             vtp_xz_idx_2 * CB_idx_2 * CM_B_idx_2;
  xB_idx_3 = (xB_idx_3 * vts_xz_idx_2 - yB_idx_3 * vr_ycomp[3]) -
             vr_xy_idx_0 * CB_idx_3 * d13;
  *Mx =
      ((((((CMx + CFx) + vthr_x) + CFy) + ThrRPMOut) + AilDefOut) + MX_thr) +
      (0.0 * F_Thr_idx_2 - -0.0058899999999999994 * F_Thr_idx_1);
  *My = (((((((ElevDefOut + vtp_xz_idx_0) +
                          ((CL_ts_idx_0 + CL_ts_idx_1) + vts_xz_idx_1)) +
                         ((zt_idx_0 + zt_idx_1) + zt_idx_2)) +
                        (((zr_idx_0 + zr_idx_1) + zr_idx_2) + zr_idx_3)) +
                       (((zB_idx_0 + zB_idx_1) + zB_idx_2) + zB_idx_3)) +
                      MX_thr_tmp * (-CMy * b_CFY_tmp + CMz * CFY_tmp) * 0.5) +
                     W_Thr_SI * 5.0E-5 * -w_r) +
                    (-0.0058899999999999994 * F_Thr_idx_0 - 0.27 * F_Thr_idx_2);
  *Mz =
      (((((((vthr_z + psiT) + RudDefOut) + ((xt_idx_0 + xt_idx_1) + xt_idx_2)) +
          (((xr_idx_0 + xr_idx_1) + xr_idx_2) + xr_idx_3)) +
         (((xB_idx_0 + xB_idx_1) + xB_idx_2) + xB_idx_3)) +
        MX_thr_tmp * (-CMy * CFY_tmp - CMz * b_CFY_tmp) * 0.5) +
       W_Thr_SI * 5.0E-5 * w_q) +
      (0.27 * F_Thr_idx_1 - 0.0 * F_Thr_idx_0);

      /* ---- pointers to output ---- */
      // *Fx = F_total_body[0];
      // *Fy = F_total_body[1];
      // *Fz = F_total_body[2];
      // *Mx = M_total_body[0];
      // *My = M_total_body[1];
      // *Mz = M_total_body[2];
}

//
// File trailer for McFoamy_FM_v3.cpp
//
// [EOF]
//
