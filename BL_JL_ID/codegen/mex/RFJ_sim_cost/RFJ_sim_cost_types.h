/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * RFJ_sim_cost_types.h
 *
 * Code generation for function 'RFJ_sim_cost'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"

/* Type Definitions */
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T
struct emxArray_real_T {
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};
#endif /* struct_emxArray_real_T */
#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T
typedef struct emxArray_real_T emxArray_real_T;
#endif /* typedef_emxArray_real_T */

#ifndef typedef_b_RFJ_sim_cost
#define typedef_b_RFJ_sim_cost
typedef struct {
  real_T zhat[26404];
} b_RFJ_sim_cost;
#endif /* typedef_b_RFJ_sim_cost */

#ifndef typedef_RFJ_sim_costStackData
#define typedef_RFJ_sim_costStackData
typedef struct {
  b_RFJ_sim_cost f0;
} RFJ_sim_costStackData;
#endif /* typedef_RFJ_sim_costStackData */

/* End of code generation (RFJ_sim_cost_types.h) */
