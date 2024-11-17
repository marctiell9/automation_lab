/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ode45.c
 *
 * Code generation for function 'ode45'
 *
 */

/* Include files */
#include "ode45.h"
#include "RFJ_sim_cost_emxutil.h"
#include "RFJ_sim_cost_mexutil.h"
#include "RFJ_sim_cost_types.h"
#include "rt_nonfinite.h"
#include "warning.h"
#include "mwmathutil.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo c_emlrtRSI = {
    17,      /* lineNo */
    "ode45", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023b\\toolbox\\eml\\lib\\matlab\\funfun\\ode45.m" /* pathName
                                                                         */
};

static emlrtRSInfo d_emlrtRSI = {
    631,                  /* lineNo */
    "explicitRungeKutta", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pathName */
};

static emlrtRSInfo e_emlrtRSI = {
    429,                  /* lineNo */
    "explicitRungeKutta", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pathName */
};

static emlrtRSInfo f_emlrtRSI = {
    416,                  /* lineNo */
    "explicitRungeKutta", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pathName */
};

static emlrtRSInfo g_emlrtRSI = {
    419,                  /* lineNo */
    "explicitRungeKutta", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pathName */
};

static emlrtRSInfo h_emlrtRSI = {
    418,                  /* lineNo */
    "explicitRungeKutta", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pathName */
};

static emlrtRSInfo i_emlrtRSI = {
    44,       /* lineNo */
    "mpower", /* fcnName */
    "C:\\Program "
    "Files\\MATLAB\\R2023b\\toolbox\\eml\\lib\\matlab\\matfun\\mpower.m" /* pathName
                                                                          */
};

static emlrtRSInfo j_emlrtRSI =
    {
        71,      /* lineNo */
        "power", /* fcnName */
        "C:\\Program "
        "Files\\MATLAB\\R2023b\\toolbox\\eml\\lib\\matlab\\ops\\power.m" /* pathName
                                                                          */
};

static emlrtRSInfo
    k_emlrtRSI =
        {
            15,        /* lineNo */
            "num2str", /* fcnName */
            "C:\\Program "
            "Files\\MATLAB\\R2023b\\toolbox\\eml\\eml\\+coder\\+"
            "internal\\num2str.m" /* pathName */
};

static emlrtMCInfo
    emlrtMCI =
        {
            53,        /* lineNo */
            19,        /* colNo */
            "flt2str", /* fName */
            "C:\\Program "
            "Files\\MATLAB\\R2023b\\toolbox\\eml\\eml\\+coder\\+"
            "internal\\flt2str.m" /* pName */
};

static emlrtRTEInfo emlrtRTEI = {
    54,                   /* lineNo */
    1,                    /* colNo */
    "explicitRungeKutta", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pName */
};

static emlrtRTEInfo b_emlrtRTEI = {
    56,                   /* lineNo */
    15,                   /* colNo */
    "explicitRungeKutta", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pName */
};

static emlrtRTEInfo c_emlrtRTEI = {
    295,                  /* lineNo */
    15,                   /* colNo */
    "explicitRungeKutta", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pName */
};

static emlrtRTEInfo e_emlrtRTEI = {
    241,                  /* lineNo */
    5,                    /* colNo */
    "explicitRungeKutta", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pName */
};

static emlrtRTEInfo f_emlrtRTEI = {
    242,                  /* lineNo */
    5,                    /* colNo */
    "explicitRungeKutta", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pName */
};

static emlrtRTEInfo g_emlrtRTEI = {
    656,                  /* lineNo */
    1,                    /* colNo */
    "explicitRungeKutta", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pName */
};

static emlrtRTEInfo h_emlrtRTEI = {
    657,                  /* lineNo */
    1,                    /* colNo */
    "explicitRungeKutta", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\explici"
    "tRungeKutta.m" /* pName */
};

static emlrtRTEInfo i_emlrtRTEI = {
    14,                  /* lineNo */
    9,                   /* colNo */
    "appendZeroColumns", /* fName */
    "C:\\Program "
    "Files\\MATLAB\\R2023b\\toolbox\\eml\\lib\\matlab\\funfun\\private\\appendZ"
    "eroColumns.m" /* pName */
};

static emlrtRSInfo
    l_emlrtRSI =
        {
            53,        /* lineNo */
            "flt2str", /* fcnName */
            "C:\\Program "
            "Files\\MATLAB\\R2023b\\toolbox\\eml\\eml\\+coder\\+"
            "internal\\flt2str.m" /* pathName */
};

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               char_T y[23]);

static const mxArray *b_sprintf(const emlrtStack *sp, const mxArray *m1,
                                const mxArray *m2, emlrtMCInfo *location);

static void emlrt_marshallIn(const emlrtStack *sp,
                             const mxArray *a__output_of_sprintf_,
                             const char_T *identifier, char_T y[23]);

static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId, char_T ret[23]);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, char_T y[23])
{
  i_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *b_sprintf(const emlrtStack *sp, const mxArray *m1,
                                const mxArray *m2, emlrtMCInfo *location)
{
  const mxArray *pArrays[2];
  const mxArray *m;
  pArrays[0] = m1;
  pArrays[1] = m2;
  return emlrtCallMATLABR2012b((emlrtConstCTX)sp, 1, &m, 2, &pArrays[0],
                               "sprintf", true, location);
}

static void emlrt_marshallIn(const emlrtStack *sp,
                             const mxArray *a__output_of_sprintf_,
                             const char_T *identifier, char_T y[23])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(a__output_of_sprintf_), &thisId, y);
  emlrtDestroyArray(&a__output_of_sprintf_);
}

static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId, char_T ret[23])
{
  static const int32_T dims[2] = {1, 23};
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "char", false, 2U,
                          (const void *)&dims[0]);
  emlrtImportCharArrayR2015b((emlrtConstCTX)sp, src, &ret[0], 23);
  emlrtDestroyArray(&src);
}

void ode45(const emlrtStack *sp, const real_T ode_workspace_th[3],
           const real_T tspan[2], const real_T b_y0[2],
           emxArray_real_T *varargout_1, emxArray_real_T *varargout_2)
{
  static const real_T x[21] = {0.2,
                               0.075,
                               0.225,
                               0.97777777777777775,
                               -3.7333333333333334,
                               3.5555555555555554,
                               2.9525986892242035,
                               -11.595793324188385,
                               9.8228928516994358,
                               -0.29080932784636487,
                               2.8462752525252526,
                               -10.757575757575758,
                               8.9064227177434727,
                               0.27840909090909088,
                               -0.2735313036020583,
                               0.091145833333333329,
                               0.0,
                               0.44923629829290207,
                               0.65104166666666663,
                               -0.322376179245283,
                               0.13095238095238096};
  static const real_T b[7] = {0.0012326388888888888,
                              0.0,
                              -0.0042527702905061394,
                              0.036979166666666667,
                              -0.05086379716981132,
                              0.0419047619047619,
                              -0.025};
  static const real_T b_b[7] = {-2.859375,
                                0.0,
                                4.0431266846361185,
                                -3.90625,
                                2.7939268867924527,
                                -1.5714285714285714,
                                1.5};
  static const real_T c_b[7] = {3.0833333333333335,
                                0.0,
                                -6.2893081761006293,
                                10.416666666666666,
                                -6.8773584905660377,
                                3.6666666666666665,
                                -4.0};
  static const real_T d_b[7] = {-1.1328125,
                                0.0,
                                2.6954177897574123,
                                -5.859375,
                                3.7610554245283021,
                                -1.9642857142857142,
                                2.5};
  static const int32_T iv[2] = {1, 7};
  static const int32_T iv1[2] = {1, 7};
  static const char_T rfmt[7] = {'%', '2', '3', '.', '1', '5', 'e'};
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack st;
  emxArray_real_T *tout;
  emxArray_real_T *yout;
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *m;
  real_T f[14];
  real_T youtnew[8];
  real_T f0[2];
  real_T fhBI2[2];
  real_T fhBI4[2];
  real_T y[2];
  real_T ystage[2];
  real_T absh;
  real_T absx;
  real_T d;
  real_T d1;
  real_T d2;
  real_T hmax;
  real_T rh;
  real_T t;
  real_T tdir;
  real_T tfinal;
  real_T *tout_data;
  real_T *varargout_1_data;
  real_T *yout_data;
  int32_T exponent;
  int32_T i;
  int32_T i1;
  int32_T j;
  int32_T k;
  int32_T ncols;
  int32_T nout;
  int32_T outidx;
  boolean_T Done;
  boolean_T MinStepExit;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  st.site = &c_emlrtRSI;
  tfinal = tspan[1];
  if (tspan[1] == 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &emlrtRTEI, "MATLAB:odearguments:TspanEndpointsNotDistinct",
        "MATLAB:odearguments:TspanEndpointsNotDistinct", 0);
  }
  MinStepExit = true;
  if (!(tspan[1] > 0.0)) {
    MinStepExit = (tspan[1] < 0.0);
  }
  if (!MinStepExit) {
    emlrtErrorWithMessageIdR2018a(&st, &b_emlrtRTEI,
                                  "MATLAB:odearguments:TspanNotMonotonic",
                                  "MATLAB:odearguments:TspanNotMonotonic", 0);
  }
  /*  Read parameters, states and inputs */
  /*  Parameters */
  /* link inertia      */
  /* for the moment it was supposed eqaul to zero, so it does not appear in the
   * model */
  /*  Torsional stiffness */
  /*  Model equations */
  /* B=zeros(1,2); */
  d = -ode_workspace_th[2] / ode_workspace_th[0];
  d1 = -ode_workspace_th[1] / ode_workspace_th[0];
  f0[0] = 0.0 * b_y0[0] + b_y0[1];
  f0[1] = b_y0[0] * d + b_y0[1] * d1;
  emxInit_real_T(&st, &tout, 2, &e_emlrtRTEI);
  i = tout->size[0] * tout->size[1];
  tout->size[0] = 1;
  tout->size[1] = 200;
  emxEnsureCapacity_real_T(&st, tout, i, &e_emlrtRTEI);
  tout_data = tout->data;
  for (i = 0; i < 200; i++) {
    tout_data[i] = 0.0;
  }
  emxInit_real_T(&st, &yout, 2, &f_emlrtRTEI);
  i = yout->size[0] * yout->size[1];
  yout->size[0] = 2;
  yout->size[1] = 200;
  emxEnsureCapacity_real_T(&st, yout, i, &f_emlrtRTEI);
  yout_data = yout->data;
  for (i = 0; i < 400; i++) {
    yout_data[i] = 0.0;
  }
  nout = 0;
  tout_data[0] = 0.0;
  yout_data[0] = b_y0[0];
  yout_data[1] = b_y0[1];
  rh = muDoubleScalarAbs(tspan[1]);
  hmax = muDoubleScalarMin(
      rh, muDoubleScalarMax(0.1 * rh, 3.5527136788005009E-15 *
                                          muDoubleScalarMax(0.0, rh)));
  if (!(hmax > 0.0)) {
    emlrtErrorWithMessageIdR2018a(&st, &c_emlrtRTEI,
                                  "MATLAB:odearguments:MaxStepLEzero",
                                  "MATLAB:odearguments:MaxStepLEzero", 0);
  }
  absh = muDoubleScalarMin(hmax, rh);
  rh = 0.0;
  for (k = 0; k < 2; k++) {
    d2 = f0[k] / muDoubleScalarMax(muDoubleScalarAbs(b_y0[k]), 0.001);
    absx = muDoubleScalarAbs(d2);
    if (muDoubleScalarIsNaN(absx) || (absx > rh)) {
      rh = absx;
    }
  }
  rh /= 0.20095091452076641;
  if (absh * rh > 1.0) {
    absh = 1.0 / rh;
  }
  absh = muDoubleScalarMax(absh, 7.90505033345994E-323);
  t = 0.0;
  y[0] = b_y0[0];
  y[1] = b_y0[1];
  memset(&f[0], 0, 14U * sizeof(real_T));
  f[0] = f0[0];
  f[1] = f0[1];
  tdir = muDoubleScalarSign(tspan[1]);
  MinStepExit = false;
  Done = false;
  int32_T exitg1;
  do {
    real_T toutnew[4];
    real_T ynew[2];
    real_T b_d2;
    real_T d3;
    real_T err;
    real_T h_tmp;
    real_T hmin;
    real_T mxerr;
    real_T tnew;
    int32_T Bcolidx;
    boolean_T NoFailedAttempts;
    exitg1 = 0;
    absx = muDoubleScalarAbs(t);
    if (muDoubleScalarIsInf(absx) || muDoubleScalarIsNaN(absx)) {
      rh = rtNaN;
    } else if (absx < 4.4501477170144028E-308) {
      rh = 4.94065645841247E-324;
    } else {
      frexp(absx, &exponent);
      rh = ldexp(1.0, exponent - 53);
    }
    hmin = 16.0 * rh;
    absh = muDoubleScalarMin(hmax, muDoubleScalarMax(hmin, absh));
    absx = tdir * absh;
    d2 = tfinal - t;
    d3 = muDoubleScalarAbs(d2);
    if (1.1 * absh >= d3) {
      absx = d2;
      absh = d3;
      Done = true;
    }
    NoFailedAttempts = true;
    int32_T exitg2;
    do {
      exitg2 = 0;
      Bcolidx = 0;
      for (j = 0; j < 5; j++) {
        Bcolidx += j;
        ystage[0] = y[0];
        ystage[1] = y[1];
        if (!(absx == 0.0)) {
          i = (j << 1) + 1;
          for (outidx = 1; outidx <= i; outidx += 2) {
            rh = absx * x[Bcolidx + ((outidx - 1) >> 1)];
            i1 = outidx + 1;
            for (k = outidx; k <= i1; k++) {
              ncols = k - outidx;
              ystage[ncols] += f[k - 1] * rh;
            }
          }
        }
        /*  Read parameters, states and inputs */
        /*  Parameters */
        /* link inertia      */
        /* for the moment it was supposed eqaul to zero, so it does not appear
         * in the model */
        /*  Torsional stiffness */
        /*  Model equations */
        /* B=zeros(1,2); */
        d2 = 0.0 * ystage[0];
        d2 += ystage[1];
        ncols = (j + 1) << 1;
        f[ncols] = d2;
        d2 = ystage[0] * d;
        d2 += ystage[1] * d1;
        f[ncols + 1] = d2;
      }
      tnew = t + absx;
      ynew[0] = y[0];
      ynew[1] = y[1];
      if (!(absx == 0.0)) {
        for (outidx = 0; outidx <= 10; outidx += 2) {
          rh = absx * x[(Bcolidx + (outidx >> 1)) + 5];
          i = outidx + 2;
          for (k = outidx + 1; k <= i; k++) {
            ncols = (k - outidx) - 1;
            ynew[ncols] += f[k - 1] * rh;
          }
        }
      }
      /*  Read parameters, states and inputs */
      /*  Parameters */
      /* link inertia      */
      /* for the moment it was supposed eqaul to zero, so it does not appear in
       * the model */
      /*  Torsional stiffness */
      /*  Model equations */
      /* B=zeros(1,2); */
      toutnew[0] = 0.0;
      toutnew[2] = 1.0;
      toutnew[1] = d;
      toutnew[3] = d1;
      for (i = 0; i < 2; i++) {
        f[i + 12] = toutnew[i] * ynew[0] + toutnew[i + 2] * ynew[1];
        d2 = 0.0;
        for (i1 = 0; i1 < 7; i1++) {
          d2 += f[i + (i1 << 1)] * b[i1];
        }
        ystage[i] = d2;
      }
      if (Done) {
        tnew = tfinal;
      }
      h_tmp = tnew - t;
      mxerr = 0.0;
      for (k = 0; k < 2; k++) {
        rh = muDoubleScalarAbs(ystage[k]);
        absx = muDoubleScalarAbs(y[k]);
        b_d2 = muDoubleScalarAbs(ynew[k]);
        if ((absx > b_d2) || muDoubleScalarIsNaN(b_d2)) {
          if (absx > 0.001) {
            rh /= absx;
          } else {
            rh /= 0.001;
          }
        } else if (b_d2 > 0.001) {
          rh /= b_d2;
        } else {
          rh /= 0.001;
        }
        if ((rh > mxerr) || muDoubleScalarIsNaN(rh)) {
          mxerr = rh;
        }
      }
      err = absh * mxerr;
      if (!(err <= 0.001)) {
        if (absh <= hmin) {
          char_T b_str[23];
          char_T str[23];
          b_st.site = &h_emlrtRSI;
          c_st.site = &k_emlrtRSI;
          b_y = NULL;
          m = emlrtCreateCharArray(2, &iv[0]);
          emlrtInitCharArrayR2013a(&c_st, 7, m, &rfmt[0]);
          emlrtAssign(&b_y, m);
          d_st.site = &l_emlrtRSI;
          emlrt_marshallIn(
              &d_st, b_sprintf(&d_st, b_y, emlrt_marshallOut(t), &emlrtMCI),
              "<output of sprintf>", str);
          b_st.site = &g_emlrtRSI;
          c_st.site = &k_emlrtRSI;
          c_y = NULL;
          m = emlrtCreateCharArray(2, &iv1[0]);
          emlrtInitCharArrayR2013a(&c_st, 7, m, &rfmt[0]);
          emlrtAssign(&c_y, m);
          d_st.site = &l_emlrtRSI;
          emlrt_marshallIn(
              &d_st, b_sprintf(&d_st, c_y, emlrt_marshallOut(hmin), &emlrtMCI),
              "<output of sprintf>", b_str);
          b_st.site = &f_emlrtRSI;
          warning(&b_st, str, b_str);
          MinStepExit = true;
          exitg2 = 1;
        } else {
          if (NoFailedAttempts) {
            NoFailedAttempts = false;
            b_st.site = &e_emlrtRSI;
            c_st.site = &i_emlrtRSI;
            d_st.site = &j_emlrtRSI;
            absh = muDoubleScalarMax(
                hmin,
                absh * muDoubleScalarMax(
                           0.1, 0.8 * muDoubleScalarPower(0.001 / err, 0.2)));
          } else {
            absh = muDoubleScalarMax(hmin, 0.5 * absh);
          }
          absx = tdir * absh;
          Done = false;
        }
      } else {
        exitg2 = 1;
      }
    } while (exitg2 == 0);
    if (MinStepExit) {
      exitg1 = 1;
    } else {
      outidx = nout + 1;
      d2 = t + h_tmp * 0.25;
      rh = d2;
      toutnew[0] = d2;
      d2 = t + h_tmp * 0.5;
      absx = d2;
      toutnew[1] = d2;
      d2 = t + h_tmp * 0.75;
      toutnew[2] = d2;
      toutnew[3] = tnew;
      for (i = 0; i < 2; i++) {
        ystage[i] = f[i] * h_tmp;
        d3 = 0.0;
        b_d2 = 0.0;
        mxerr = 0.0;
        for (i1 = 0; i1 < 7; i1++) {
          hmin = f[i + (i1 << 1)];
          d3 += hmin * (h_tmp * b_b[i1]);
          b_d2 += hmin * (h_tmp * c_b[i1]);
          mxerr += hmin * (h_tmp * d_b[i1]);
        }
        fhBI4[i] = mxerr;
        f0[i] = b_d2;
        fhBI2[i] = d3;
      }
      __m128d r;
      __m128d r1;
      __m128d r2;
      __m128d r3;
      __m128d r4;
      __m128d r5;
      r = _mm_loadu_pd(&fhBI4[0]);
      r1 = _mm_loadu_pd(&f0[0]);
      r2 = _mm_loadu_pd(&fhBI2[0]);
      r3 = _mm_loadu_pd(&ystage[0]);
      r4 = _mm_loadu_pd(&y[0]);
      r5 = _mm_set1_pd((rh - t) / h_tmp);
      _mm_storeu_pd(
          &youtnew[0],
          _mm_add_pd(
              _mm_mul_pd(
                  _mm_add_pd(
                      _mm_mul_pd(
                          _mm_add_pd(
                              _mm_mul_pd(_mm_add_pd(_mm_mul_pd(r, r5), r1), r5),
                              r2),
                          r5),
                      r3),
                  r5),
              r4));
      r = _mm_loadu_pd(&fhBI4[0]);
      r1 = _mm_loadu_pd(&f0[0]);
      r2 = _mm_loadu_pd(&fhBI2[0]);
      r3 = _mm_loadu_pd(&ystage[0]);
      r4 = _mm_loadu_pd(&y[0]);
      r5 = _mm_set1_pd((absx - t) / h_tmp);
      _mm_storeu_pd(
          &youtnew[2],
          _mm_add_pd(
              _mm_mul_pd(
                  _mm_add_pd(
                      _mm_mul_pd(
                          _mm_add_pd(
                              _mm_mul_pd(_mm_add_pd(_mm_mul_pd(r, r5), r1), r5),
                              r2),
                          r5),
                      r3),
                  r5),
              r4));
      r = _mm_loadu_pd(&fhBI4[0]);
      r1 = _mm_loadu_pd(&f0[0]);
      r2 = _mm_loadu_pd(&fhBI2[0]);
      r3 = _mm_loadu_pd(&ystage[0]);
      r4 = _mm_loadu_pd(&y[0]);
      r5 = _mm_set1_pd((d2 - t) / h_tmp);
      _mm_storeu_pd(
          &youtnew[4],
          _mm_add_pd(
              _mm_mul_pd(
                  _mm_add_pd(
                      _mm_mul_pd(
                          _mm_add_pd(
                              _mm_mul_pd(_mm_add_pd(_mm_mul_pd(r, r5), r1), r5),
                              r2),
                          r5),
                      r3),
                  r5),
              r4));
      youtnew[6] = ynew[0];
      youtnew[7] = ynew[1];
      nout += 4;
      if (nout + 1 > tout->size[1]) {
        ncols = tout->size[1];
        i = tout->size[0] * tout->size[1];
        tout->size[0] = 1;
        tout->size[1] += 200;
        emxEnsureCapacity_real_T(&st, tout, i, &i_emlrtRTEI);
        tout_data = tout->data;
        i = yout->size[0] * yout->size[1];
        yout->size[0] = 2;
        yout->size[1] += 200;
        emxEnsureCapacity_real_T(&st, yout, i, &i_emlrtRTEI);
        yout_data = yout->data;
        for (j = 0; j < 200; j++) {
          Bcolidx = ncols + j;
          tout_data[Bcolidx] = 0.0;
          yout_data[2 * Bcolidx] = 0.0;
          yout_data[2 * Bcolidx + 1] = 0.0;
        }
      }
      for (k = 0; k < 4; k++) {
        Bcolidx = k + outidx;
        tout_data[Bcolidx] = toutnew[k];
        ncols = k << 1;
        yout_data[2 * Bcolidx] = youtnew[ncols];
        yout_data[2 * Bcolidx + 1] = youtnew[ncols + 1];
      }
      if (Done) {
        exitg1 = 1;
      } else {
        if (NoFailedAttempts) {
          b_st.site = &d_emlrtRSI;
          c_st.site = &i_emlrtRSI;
          d_st.site = &j_emlrtRSI;
          rh = 1.25 * muDoubleScalarPower(err / 0.001, 0.2);
          if (rh > 0.2) {
            absh /= rh;
          } else {
            absh *= 5.0;
          }
        }
        t = tnew;
        y[0] = ynew[0];
        f[0] = f[12];
        y[1] = ynew[1];
        f[1] = f[13];
      }
    }
  } while (exitg1 == 0);
  if (nout + 1 < 1) {
    ncols = -1;
  } else {
    ncols = nout;
  }
  i = varargout_1->size[0];
  varargout_1->size[0] = ncols + 1;
  emxEnsureCapacity_real_T(&st, varargout_1, i, &g_emlrtRTEI);
  varargout_1_data = varargout_1->data;
  for (i = 0; i <= ncols; i++) {
    varargout_1_data[i] = tout_data[i];
  }
  emxFree_real_T(&st, &tout);
  i = varargout_2->size[0] * varargout_2->size[1];
  varargout_2->size[0] = ncols + 1;
  varargout_2->size[1] = 2;
  emxEnsureCapacity_real_T(&st, varargout_2, i, &h_emlrtRTEI);
  tout_data = varargout_2->data;
  for (i = 0; i < 2; i++) {
    for (i1 = 0; i1 <= ncols; i1++) {
      tout_data[i1 + varargout_2->size[0] * i] = yout_data[i + 2 * i1];
    }
  }
  emxFree_real_T(&st, &yout);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (ode45.c) */
