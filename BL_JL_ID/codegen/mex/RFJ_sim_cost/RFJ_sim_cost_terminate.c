/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * RFJ_sim_cost_terminate.c
 *
 * Code generation for function 'RFJ_sim_cost_terminate'
 *
 */

/* Include files */
#include "RFJ_sim_cost_terminate.h"
#include "RFJ_sim_cost_data.h"
#include "_coder_RFJ_sim_cost_mex.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void RFJ_sim_cost_atexit(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void RFJ_sim_cost_terminate(void)
{
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (RFJ_sim_cost_terminate.c) */
