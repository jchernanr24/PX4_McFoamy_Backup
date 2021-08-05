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
