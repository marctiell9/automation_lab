/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * RFJ_sim_cost_initialize.c
 *
 * Code generation for function 'RFJ_sim_cost_initialize'
 *
 */

/* Include files */
#include "RFJ_sim_cost_initialize.h"
#include "RFJ_sim_cost_data.h"
#include "_coder_RFJ_sim_cost_mex.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void RFJ_sim_cost_once(void);

/* Function Definitions */
static void RFJ_sim_cost_once(void)
{
  mex_InitInfAndNan();
}

void RFJ_sim_cost_initialize(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2022b(&st);
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  if (emlrtFirstTimeR2012b(emlrtRootTLSGlobal)) {
    RFJ_sim_cost_once();
  }
}

/* End of code generation (RFJ_sim_cost_initialize.c) */
