/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ode45.h
 *
 * Code generation for function 'ode45'
 *
 */

#pragma once

/* Include files */
#include "RFJ_sim_cost_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void ode45(const emlrtStack *sp, const real_T ode_workspace_th[3],
           const real_T tspan[2], const real_T b_y0[2],
           emxArray_real_T *varargout_1, emxArray_real_T *varargout_2);

/* End of code generation (ode45.h) */
