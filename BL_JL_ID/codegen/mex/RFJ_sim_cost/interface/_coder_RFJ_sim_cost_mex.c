/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_RFJ_sim_cost_mex.c
 *
 * Code generation for function '_coder_RFJ_sim_cost_mex'
 *
 */

/* Include files */
#include "_coder_RFJ_sim_cost_mex.h"
#include "RFJ_sim_cost_data.h"
#include "RFJ_sim_cost_initialize.h"
#include "RFJ_sim_cost_terminate.h"
#include "RFJ_sim_cost_types.h"
#include "_coder_RFJ_sim_cost_api.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void RFJ_sim_cost_mexFunction(RFJ_sim_costStackData *SD, int32_T nlhs,
                              mxArray *plhs[2], int32_T nrhs,
                              const mxArray *prhs[8])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs[2];
  int32_T i;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 8) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 8, 4,
                        12, "RFJ_sim_cost");
  }
  if (nlhs > 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 12,
                        "RFJ_sim_cost");
  }
  /* Call the function. */
  RFJ_sim_cost_api(SD, prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    i = 1;
  } else {
    i = nlhs;
  }
  emlrtReturnArrays(i, &plhs[0], &outputs[0]);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  RFJ_sim_costStackData *RFJ_sim_costStackDataGlobal = NULL;
  RFJ_sim_costStackDataGlobal = (RFJ_sim_costStackData *)emlrtMxCalloc(
      (size_t)1, (size_t)1U * sizeof(RFJ_sim_costStackData));
  mexAtExit(&RFJ_sim_cost_atexit);
  /* Module initialization. */
  RFJ_sim_cost_initialize();
  /* Dispatch the entry-point. */
  RFJ_sim_cost_mexFunction(RFJ_sim_costStackDataGlobal, nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  RFJ_sim_cost_terminate();
  emlrtMxFree(RFJ_sim_costStackDataGlobal);
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "windows-1252", true);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_RFJ_sim_cost_mex.c) */
