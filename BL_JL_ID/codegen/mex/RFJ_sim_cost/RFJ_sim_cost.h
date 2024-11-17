/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * RFJ_sim_cost.h
 *
 * Code generation for function 'RFJ_sim_cost'
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
void RFJ_sim_cost(RFJ_sim_costStackData *SD, const emlrtStack *sp,
                  const real_T x[2], const real_T z0[2],
                  const real_T uin[13202], const real_T ymeas[13202], real_T th,
                  real_T Ts, real_T Q, const real_T scaling[2], real_T *cost,
                  real_T ysim[13202]);

/* End of code generation (RFJ_sim_cost.h) */
