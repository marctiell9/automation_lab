/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_RFJ_sim_cost_info.c
 *
 * Code generation for function 'RFJ_sim_cost'
 *
 */

/* Include files */
#include "_coder_RFJ_sim_cost_info.h"
#include "emlrt.h"
#include "tmwtypes.h"

/* Function Declarations */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void);

/* Function Definitions */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void)
{
  const mxArray *nameCaptureInfo;
  const char_T *data[6] = {
      "789ced55cd4edb4010fe5c68858468d34bafe50950851451f5d6a44d6b94849f18a5a241"
      "c18d97b2e0b523ff44c04bf439fa0a7d04c4812b9c382271e405183b"
      "ebfc585dc508ea2a55578a67befd66f71b4f461e687a4d03f01cfd75bed0b7d2a020ed13"
      "8caf34afa5e2b4f1703cc5ecd8b984ff216dc775027614f481630a36",
      "3869b9823ba61318c75d068ff9aedd6356ccec719b195cb0c628a847485446a80188a8c8"
      "2fefb3ce612314f0f6fd6186f62818d4e38df6fbf79dcd588f3d453d"
      "5e4a7b26ed577cc40eca788716b6e0835e959e2d0898e475e092ff81767d1c2220d4251c"
      "c201474fc672dabfa2dd22fda2f3228e0e886358a6bdf774af8135d4",
      "c833a09357471b55422562755814e7503ca79c39299ab1efd25e8b22aa14bb1a3f75caa3"
      "854d5408b7635d4136ca30525b22345abfaea23e59eb37afa85f21c5"
      "6f5656db3e176de67963fabb0fd47fa6d4ef33961b7eb3d9e3f50b4be16425fd722aedb4"
      "f74b948337e896c9ffd78b8cf54bdb61fc5c6c6f8e2fb53cf56e2f8a",
      "7e9e7ac9fa5b7a478afbb2f6ff2b855e21c5d7c4ca4ab97ad0ac7bcdb0d1fbce6da3b3fd"
      "b634cc637d82cea43ca0c079ddffff3bf267be2359eb37a3a85f4132"
      "346f1e55efbe7366f181fdb1a3d04bfae397b4d3d91ff9cf95d761be7365fbfaa79ea75e"
      "b2fef5b9d2d868ba5d4b1c7c5e334f4a453730beac1bfaa7e99f2b77",
      "0ddbf099",
      ""};
  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(&data[0], 3560U, &nameCaptureInfo);
  return nameCaptureInfo;
}

mxArray *emlrtMexFcnProperties(void)
{
  mxArray *xEntryPoints;
  mxArray *xInputs;
  mxArray *xResult;
  const char_T *propFieldName[9] = {"Version",
                                    "ResolvedFunctions",
                                    "Checksum",
                                    "EntryPoints",
                                    "CoverageInfo",
                                    "IsPolymorphic",
                                    "PropertyList",
                                    "UUID",
                                    "ClassEntryPointIsHandle"};
  const char_T *epFieldName[8] = {
      "Name",     "NumberOfInputs", "NumberOfOutputs", "ConstantInputs",
      "FullPath", "TimeStamp",      "Constructor",     "Visible"};
  xEntryPoints =
      emlrtCreateStructMatrix(1, 1, 8, (const char_T **)&epFieldName[0]);
  xInputs = emlrtCreateLogicalMatrix(1, 8);
  emlrtSetField(xEntryPoints, 0, "Name", emlrtMxCreateString("RFJ_sim_cost"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(8.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(2.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 0, "FullPath",
      emlrtMxCreateString(
          "C:"
          "\\Users\\marco\\Desktop\\universit\xc3\xa0\\5\\semestre2\\AUTOMATION"
          "_LAB\\Identification\\BL_JL_ID\\RFJ_sim_cost.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(739312.70571759262));
  emlrtSetField(xEntryPoints, 0, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 0, "Visible", emlrtMxCreateLogicalScalar(true));
  xResult =
      emlrtCreateStructMatrix(1, 1, 9, (const char_T **)&propFieldName[0]);
  emlrtSetField(xResult, 0, "Version",
                emlrtMxCreateString("23.2.0.2391609 (R2023b) Update 2"));
  emlrtSetField(xResult, 0, "ResolvedFunctions",
                (mxArray *)c_emlrtMexFcnResolvedFunctionsI());
  emlrtSetField(xResult, 0, "Checksum",
                emlrtMxCreateString("La13huTwz6WtTdJOtEtTmC"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/* End of code generation (_coder_RFJ_sim_cost_info.c) */
