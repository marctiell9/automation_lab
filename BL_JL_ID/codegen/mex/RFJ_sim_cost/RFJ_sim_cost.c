/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * RFJ_sim_cost.c
 *
 * Code generation for function 'RFJ_sim_cost'
 *
 */

/* Include files */
#include "RFJ_sim_cost.h"
#include "RFJ_sim_cost_data.h"
#include "RFJ_sim_cost_emxutil.h"
#include "RFJ_sim_cost_types.h"
#include "ode45.h"
#include "rt_nonfinite.h"
#include "sumMatrixIncludeNaN.h"
#include <emmintrin.h>
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = {
    9,              /* lineNo */
    "RFJ_sim_cost", /* fcnName */
    "C:\\Users\\marco\\Desktop\\universit\xc3\xa0\\5\\semestre2\\AUTOMATION_"
    "LAB\\Identification\\BL_JL_ID\\RFJ_sim_cost.m" /* pathName */
};

static emlrtRSInfo b_emlrtRSI = {
    21,            /* lineNo */
    "RFJ_sim_err", /* fcnName */
    "C:\\Users\\marco\\Desktop\\universit\xc3\xa0\\5\\semestre2\\AUTOMATION_"
    "LAB\\Identification\\BL_JL_ID\\RFJ_sim_err.m" /* pathName */
};

static emlrtBCInfo emlrtBCI = {
    -1,            /* iFirst */
    -1,            /* iLast */
    22,            /* lineNo */
    42,            /* colNo */
    "zhat_temp",   /* aName */
    "RFJ_sim_err", /* fName */
    "C:\\Users\\marco\\Desktop\\universit\xc3\xa0\\5\\semestre2\\AUTOMATION_"
    "LAB\\Identification\\BL_JL_ID\\RFJ_sim_err.m", /* pName */
    0                                               /* checkKind */
};

static emlrtRTEInfo d_emlrtRTEI = {
    1,              /* lineNo */
    32,             /* colNo */
    "RFJ_sim_cost", /* fName */
    "C:\\Users\\marco\\Desktop\\universit\xc3\xa0\\5\\semestre2\\AUTOMATION_"
    "LAB\\Identification\\BL_JL_ID\\RFJ_sim_cost.m" /* pName */
};

/* Function Definitions */
void RFJ_sim_cost(RFJ_sim_costStackData *SD, const emlrtStack *sp,
                  const real_T x[2], const real_T z0[2],
                  const real_T uin[13202], const real_T ymeas[13202], real_T th,
                  real_T Ts, real_T Q, const real_T scaling[2], real_T *cost,
                  real_T ysim[13202])
{
  __m128d r;
  __m128d r1;
  emlrtStack b_st;
  emlrtStack st;
  emxArray_real_T *t;
  emxArray_real_T *zhat_temp;
  real_T err_vec[13202];
  real_T b_th[3];
  real_T b_x[2];
  real_T dv[2];
  real_T s;
  real_T *zhat_temp_data;
  int32_T ind;
  int32_T zhat_tmp;
  (void)uin;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  /*  Function that computes the output trajectory of the vehicle model */
  /*  introduced in Lab. session A and the quadratic cost given by */
  /*  the sum of squared weighted errors between the simulated output */
  /*  and the measured one. Variable x corresponds to the vehicle parameters  */
  /*  to be estimated, possibly scaled. The simulation is carried out with
   * forward finite */
  /*  differences */
  st.site = &emlrtRSI;
  /*  Model parameters */
  r = _mm_loadu_pd(&x[0]);
  r1 = _mm_loadu_pd(&scaling[0]);
  _mm_storeu_pd(&b_x[0], _mm_div_pd(r, r1));
  b_th[0] = th;
  b_th[1] = b_x[0];
  b_th[2] = b_x[1];
  /*  Initialize simulation output */
  /*  number of samples */
  memset(&SD->f0.zhat[0], 0, 26404U * sizeof(real_T));
  /*  matrix with states */
  SD->f0.zhat[0] = z0[0];
  SD->f0.zhat[1] = z0[1];
  /*  FFD */
  /*  for ind=2:N */
  /*      zdot               =   RFJ(0,zhat(:,ind-1),uin(:,ind-1),th); */
  /*      zhat(:,ind)    =   zhat(:,ind-1)+Ts*zdot; */
  /*  end */
  /*  ODE45 */
  dv[0] = 0.0;
  dv[1] = Ts;
  emxInit_real_T(&st, &t, 1, &d_emlrtRTEI);
  emxInit_real_T(&st, &zhat_temp, 2, &d_emlrtRTEI);
  for (ind = 0; ind < 13201; ind++) {
    b_st.site = &b_emlrtRSI;
    ode45(&b_st, b_th, dv, &SD->f0.zhat[ind << 1], t, zhat_temp);
    zhat_temp_data = zhat_temp->data;
    if (zhat_temp->size[0] < 1) {
      emlrtDynamicBoundsCheckR2012b(zhat_temp->size[0], 1, zhat_temp->size[0],
                                    &emlrtBCI, &st);
    }
    zhat_tmp = (ind + 1) << 1;
    SD->f0.zhat[zhat_tmp] = zhat_temp_data[zhat_temp->size[0] - 1];
    SD->f0.zhat[zhat_tmp + 1] =
        zhat_temp_data[(zhat_temp->size[0] + zhat_temp->size[0]) - 1];
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(&st);
    }
  }
  emxFree_real_T(&st, &zhat_temp);
  emxFree_real_T(&st, &t);
  /*  Collect simulated output  */
  /*  Compute weighted errors */
  /*  Stack errors in one vector */
  for (zhat_tmp = 0; zhat_tmp < 13202; zhat_tmp++) {
    s = SD->f0.zhat[zhat_tmp << 1];
    ysim[zhat_tmp] = s;
    err_vec[zhat_tmp] = Q * (ymeas[zhat_tmp] - s) / 114.8999564838908;
  }
  /*  Compute sum of squared errors */
  for (zhat_tmp = 0; zhat_tmp <= 13200; zhat_tmp += 2) {
    r = _mm_loadu_pd(&err_vec[zhat_tmp]);
    _mm_storeu_pd(&err_vec[zhat_tmp], _mm_mul_pd(r, r));
  }
  s = sumColumnB4(err_vec);
  s += b_sumColumnB4(err_vec, 4097);
  s += b_sumColumnB4(err_vec, 8193);
  *cost = s + sumColumnB(err_vec);
  /* questa è la f(x) in tale problema */
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (RFJ_sim_cost.c) */
