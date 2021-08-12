//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: main.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 05-Aug-2021 13:41:18
//

/*************************************************************************/
/* This automatically generated example C++ main file shows how to call  */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/

// Include Files
#include "main.h"
#include "McFoamy_FM_v3.h"
#include "McFoamy_FM_v3_terminate.h"
#include "rt_nonfinite.h"

// Function Declarations
static double argInit_real_T();

static void main_McFoamy_FM_v3();

// Function Definitions
//
// Arguments    : void
// Return Type  : double
//
static double argInit_real_T()
{
  return 0.0;
}

//
// Arguments    : void
// Return Type  : void
//
static void main_McFoamy_FM_v3()
{
  double AER_lim[3];
  double F_total_body[3];
  double M_total_body[3];
  double Ail_def_tmp;
  // Initialize function 'McFoamy_FM_v3' input arguments.
  Ail_def_tmp = argInit_real_T();
  // Call the entry-point 'McFoamy_FM_v3'.
  McFoamy_FM_v3(Ail_def_tmp, Ail_def_tmp, Ail_def_tmp, Ail_def_tmp, Ail_def_tmp,
                Ail_def_tmp, Ail_def_tmp, Ail_def_tmp, Ail_def_tmp, Ail_def_tmp,
                F_total_body, M_total_body, AER_lim);
}

//
// Arguments    : int argc
//                char **argv
// Return Type  : int
//
int main(int, char **)
{
  // The initialize function is being called automatically from your entry-point
  // function. So, a call to initialize is not included here. Invoke the
  // entry-point functions.
  // You can call entry-point functions multiple times.
  main_McFoamy_FM_v3();
  // Terminate the application.
  // You do not need to do this more than one time.
  McFoamy_FM_v3_terminate();
  return 0;
}

//
// File trailer for main.cpp
//
// [EOF]
//
