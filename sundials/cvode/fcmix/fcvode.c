/*
 * -----------------------------------------------------------------
 * $Revision: 1.61 $
 * $Date: 2006/02/10 00:03:09 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Alan C. Hindmarsh, Radu Serban and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see sundials/cvode/LICENSE.
 * -----------------------------------------------------------------
 * This is the implementation file for the Fortran interface to
 * the CVODE package.  See fcvode.h for usage.
 * NOTE: some routines are necessarily stored elsewhere to avoid
 * linking problems.  Therefore, see also fcvpreco.c, fcvpsol.c,
 * and fcvjtimes.c for all the options available.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fcvode.h"           /* actual function names, prototypes, global vars. */
#include "cvode.h"            /* CVODE constants and prototypes                  */

#include "cvode_band.h"       /* prototypes for CVBAND interface routines        */
#include "cvode_dense.h"      /* prototypes for CVDENSE interface routines       */
#include "cvode_diag.h"       /* prototypes for CVDIAG interface routines        */
#include "cvode_spgmr.h"      /* prototypes for CVSPGMR interface routines       */
#include "cvode_spbcgs.h"     /* prototypes for CVSPBCG interface routines       */
#include "cvode_sptfqmr.h"    /* prototypes for CVSPTFQMR interface routines     */

#include "cvode_impl.h"       /* definition of CVodeMem type                     */

#include "sundials_nvector.h" /* definitions of type N_Vector and vector macros  */
#include "sundials_types.h"   /* definition of type realtype                     */

/***************************************************************************/

/* Definitions for global variables shared amongst various routines */

void *CV_cvodemem;
long int *CV_iout;
realtype *CV_rout;
int CV_nrtfn;
int CV_ls;

/***************************************************************************/

/* private constant(s) */
#define ZERO RCONST(0.0)

/***************************************************************************/

/* Prototypes of the Fortran routines */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
  extern void FCV_FUN(realtype*,     /* T    */
                      realtype*,     /* Y    */
                      realtype*,     /* YDOT */
                      long int*,     /* IPAR */
                      realtype*,     /* RPAR */
                      int*);         /* IER  */
#ifdef __cplusplus
}
#endif

/**************************************************************************/

void FCV_MALLOC(realtype *t0, realtype *y0, 
                int *meth, int *itmeth, int *iatol, 
                realtype *rtol, realtype *atol,
                long int *iout, realtype *rout,
                long int *ipar, realtype *rpar,
                int *ier)
{
  int lmm, iter, itol;
  N_Vector Vatol;
  void *atolptr;
  FCVUserData CV_userdata;

  *ier = 0;

  /* Check for required vector operations */
  if(F2C_CVODE_vec->ops->nvgetarraypointer == NULL ||
     F2C_CVODE_vec->ops->nvsetarraypointer == NULL) {
    *ier = -1;
    printf("A required vector operation is not implemented.\n\n");
    return;
  }

  /* Initialize all pointers to NULL */
  CV_cvodemem = NULL;
  Vatol = NULL;
  atolptr = NULL;

  /* Create CVODE object */
  lmm = (*meth == 1) ? CV_ADAMS : CV_BDF;
  iter = (*itmeth == 1) ? CV_FUNCTIONAL : CV_NEWTON;
  CV_cvodemem = CVodeCreate(lmm, iter);
  if (CV_cvodemem == NULL) {
    *ier = -1;
    return;
  }

  /* Set and attach user data */
  CV_userdata = NULL;
  CV_userdata = (FCVUserData) malloc(sizeof *CV_userdata);
  if (CV_userdata == NULL) {
    *ier = -1;
    return;
  }
  CV_userdata->rpar = rpar;
  CV_userdata->ipar = ipar;

  *ier = CVodeSetFdata(CV_cvodemem, CV_userdata);
  if(*ier != CV_SUCCESS) {
    free(CV_userdata); CV_userdata = NULL;
    *ier = -1;
    return;
  }

  /* Set data in F2C_CVODE_vec to y0 */
  N_VSetArrayPointer(y0, F2C_CVODE_vec);

  /* Treat absolute tolerances */
  itol = -1;
  switch (*iatol) {
  case 1:
    itol = CV_SS; 
    atolptr = (void *) atol; 
    break;
  case 2:
    itol = CV_SV; 
    Vatol = NULL;
    Vatol = N_VCloneEmpty(F2C_CVODE_vec);
    if (Vatol == NULL) {
      free(CV_userdata); CV_userdata = NULL;
      *ier = -1;
      return;
    }
    N_VSetArrayPointer(atol, Vatol);
    atolptr = (void *) Vatol; 
    break;
  case 3:
    itol = CV_WF;
    break;
  }

  /* Call CVodeMalloc */
  *ier = CVodeMalloc(CV_cvodemem, FCVf, *t0, F2C_CVODE_vec, itol, *rtol, atolptr);

  /* destroy Vatol if allocated */
  if (itol == CV_SV) N_VDestroy(Vatol);

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_CVODE_vec);

  /* On failure, exit */
  if(*ier != CV_SUCCESS) {
    free(CV_userdata); CV_userdata = NULL;
    *ier = -1;
    return;
  }

  /* Grab optional output arrays and store them in global variables */
  CV_iout = iout;
  CV_rout = rout;

  /* Store the unit roundoff in rout for user access */
  CV_rout[5] = UNIT_ROUNDOFF;

  return;
}

/***************************************************************************/

void FCV_REINIT(realtype *t0, realtype *y0, 
                int *iatol, realtype *rtol, realtype *atol, 
                int *ier)
{
  int itol;
  N_Vector Vatol;
  void *atolptr;

  *ier = 0;

  /* Initialize all pointers to NULL */
  Vatol = NULL;
  atolptr = NULL;

  /* Set data in F2C_CVODE_vec to y0 */
  N_VSetArrayPointer(y0, F2C_CVODE_vec);

  /* Treat absolute tolerances */
  itol = -1;
  switch (*iatol) {
  case 1:
    itol = CV_SS; 
    atolptr = (void *) atol; 
    break;
  case 2:
    itol = CV_SV; 
    Vatol = NULL;
    Vatol = N_VCloneEmpty(F2C_CVODE_vec);
    if (Vatol == NULL) {
      *ier = -1;
      return;
    }
    N_VSetArrayPointer(atol, Vatol);
    atolptr = (void *) Vatol; 
    break;
  case 3:
    itol = CV_WF;
    break;
  }

  /* Call CVReInit */
  *ier = CVodeReInit(CV_cvodemem, FCVf, *t0, F2C_CVODE_vec, itol, *rtol, atolptr);

  /* destroy Vatol if allocated */
  if (itol == CV_SV) N_VDestroy(Vatol);

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_CVODE_vec);

  /* On failure, exit */
  if (*ier != CV_SUCCESS) {
    *ier = -1;
    return;
  }

  return;
}

/***************************************************************************/

void FCV_SETIIN(char key_name[], long int *ival, int *ier, int key_len)
{
  if (!strncmp(key_name,"MAX_ORD", (size_t)key_len)) 
    *ier = CVodeSetMaxOrd(CV_cvodemem, (int) *ival);
  else if (!strncmp(key_name,"MAX_NSTEPS", (size_t)key_len)) 
    *ier = CVodeSetMaxNumSteps(CV_cvodemem, (int) *ival);
  else if (!strncmp(key_name,"MAX_ERRFAIL", (size_t)key_len)) 
    *ier = CVodeSetMaxErrTestFails(CV_cvodemem, (int) *ival);
  else if (!strncmp(key_name,"MAX_NITERS", (size_t)key_len)) 
    *ier = CVodeSetMaxNonlinIters(CV_cvodemem, (int) *ival);
  else if (!strncmp(key_name,"MAX_CONVFAIL", (size_t)key_len)) 
    *ier = CVodeSetMaxConvFails(CV_cvodemem, (int) *ival);
  else if (!strncmp(key_name,"HNIL_WARNS", (size_t)key_len)) 
    *ier = CVodeSetMaxHnilWarns(CV_cvodemem, (int) *ival);
  else if (!strncmp(key_name,"STAB_LIM", (size_t)key_len)) 
    *ier = CVodeSetStabLimDet(CV_cvodemem, (int) *ival);
  else {
    *ier = -99;
    printf("FCVSETIIN: Unrecognized key.\n\n");
  }

}

/***************************************************************************/

void FCV_SETRIN(char key_name[], realtype *rval, int *ier, int key_len)
{
  if (!strncmp(key_name,"INIT_STEP", (size_t)key_len)) 
    *ier = CVodeSetInitStep(CV_cvodemem, *rval);
  else if (!strncmp(key_name,"MAX_STEP", (size_t)key_len)) 
    *ier = CVodeSetMaxStep(CV_cvodemem, *rval);
  else if (!strncmp(key_name,"MIN_STEP", (size_t)key_len)) 
    *ier = CVodeSetMinStep(CV_cvodemem, *rval);
  else if (!strncmp(key_name,"STOP_TIME", (size_t)key_len)) 
    *ier = CVodeSetStopTime(CV_cvodemem, *rval);
  else if (!strncmp(key_name,"NLCONV_COEF", (size_t)key_len)) 
    *ier = CVodeSetNonlinConvCoef(CV_cvodemem, *rval);
  else {
    *ier = -99;
    printf("FCVSETRIN: Unrecognized key.\n\n");
  }

}

/***************************************************************************/

void FCV_DENSE(long int *neq, int *ier)
{
  /* neq  is the problem size */

  *ier = CVDense(CV_cvodemem, *neq);

  CV_ls = CV_LS_DENSE;
}

/***************************************************************************/

void FCV_BAND(long int *neq, long int *mupper, long int *mlower, int *ier)
{
  /* 
     neq        is the problem size
     mupper     is the upper bandwidth
     mlower     is the lower bandwidth 
  */

  *ier = CVBand(CV_cvodemem, *neq, *mupper, *mlower);

  CV_ls = CV_LS_BAND;
}

/***************************************************************************/

void FCV_DIAG(int *ier)
{
  *ier = CVDiag(CV_cvodemem);

  CV_ls = CV_LS_DIAG;
}

/***************************************************************************/

void FCV_SPGMR(int *pretype, int *gstype, int *maxl, realtype *delt, int *ier)
{
  /* 
     pretype    the preconditioner type
     maxl       the maximum Krylov dimension
     gstype     the Gram-Schmidt process type
     delt       the linear convergence tolerance factor 
  */

  *ier = CVSpgmr(CV_cvodemem, *pretype, *maxl);
  if (*ier != CVSPILS_SUCCESS) return;

  *ier = CVSpilsSetGSType(CV_cvodemem, *gstype);
  if (*ier != CVSPILS_SUCCESS) return;

  *ier = CVSpilsSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPILS_SUCCESS) return;

  CV_ls = CV_LS_SPGMR;
}

/***************************************************************************/

void FCV_SPBCG(int *pretype, int *maxl, realtype *delt, int *ier)
{
  /* 
     pretype    the preconditioner type
     maxl       the maximum Krylov dimension
     delt       the linear convergence tolerance factor 
  */

  *ier = CVSpbcg(CV_cvodemem, *pretype, *maxl);
  if (*ier != CVSPILS_SUCCESS) return;

  *ier = CVSpilsSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPILS_SUCCESS) return;

  CV_ls = CV_LS_SPBCG;
}

/***************************************************************************/

void FCV_SPTFQMR(int *pretype, int *maxl, realtype *delt, int *ier)
{
  /* 
     pretype    the preconditioner type
     maxl       the maximum Krylov dimension
     delt       the linear convergence tolerance factor 
  */

  *ier = CVSptfqmr(CV_cvodemem, *pretype, *maxl);
  if (*ier != CVSPILS_SUCCESS) return;

  *ier = CVSpilsSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPILS_SUCCESS) return;

  CV_ls = CV_LS_SPTFQMR;
}

/***************************************************************************/

void FCV_SPGMRREINIT(int *pretype, int *gstype, realtype *delt, int *ier)
{
  /* 
     pretype    the preconditioner type
     gstype     the Gram-Schmidt process type
     delt       the linear convergence tolerance factor 
  */

  *ier = CVSpilsSetPrecType(CV_cvodemem, *pretype);
  if (*ier != CVSPILS_SUCCESS) return;

  *ier = CVSpilsSetGSType(CV_cvodemem, *gstype);
  if (*ier != CVSPILS_SUCCESS) return;

  *ier = CVSpilsSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPILS_SUCCESS) return;

  CV_ls = CV_LS_SPGMR;
}

/***************************************************************************/

void FCV_SPBCGREINIT(int *pretype, int *maxl, realtype *delt, int *ier)
{
  /* 
     pretype    the preconditioner type
     maxl       the maximum Krylov subspace dimension
     delt       the linear convergence tolerance factor 
  */

  *ier = CVSpilsSetPrecType(CV_cvodemem, *pretype);
  if (*ier != CVSPILS_SUCCESS) return;

  *ier = CVSpilsSetMaxl(CV_cvodemem, *maxl);
  if (*ier != CVSPILS_SUCCESS) return;

  *ier = CVSpilsSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPILS_SUCCESS) return;

  CV_ls = CV_LS_SPBCG;
}

/***************************************************************************/

void FCV_SPTFQMRREINIT(int *pretype, int *maxl, realtype *delt, int *ier)
{
  /* 
     pretype    the preconditioner type
     maxl       the maximum Krylov subspace dimension
     delt       the linear convergence tolerance factor 
  */

  *ier = CVSpilsSetPrecType(CV_cvodemem, *pretype);
  if (*ier != CVSPILS_SUCCESS) return;

  *ier = CVSpilsSetMaxl(CV_cvodemem, *maxl);
  if (*ier != CVSPILS_SUCCESS) return;

  *ier = CVSpilsSetDelt(CV_cvodemem, *delt);
  if (*ier != CVSPILS_SUCCESS) return;

  CV_ls = CV_LS_SPTFQMR;
}

/***************************************************************************/

void FCV_CVODE(realtype *tout, realtype *t, realtype *y, int *itask, int *ier)
{
  /* 
     tout          is the t value where output is desired
     F2C_CVODE_vec is the N_Vector containing the solution on return
     t             is the returned independent variable value
     itask         is the task indicator (1 = CV_NORMAL, 2 = CV_ONE_STEP, 
                                          3 = CV_NORMAL_TSTOP, 4 = CV_ONE_STEP_TSTOP) 
  */

  N_VSetArrayPointer(y, F2C_CVODE_vec);

  *ier = CVode(CV_cvodemem, *tout, F2C_CVODE_vec, t, *itask);

  N_VSetArrayPointer(NULL, F2C_CVODE_vec);

  /* Load optional outputs in iout & rout */
  CVodeGetWorkSpace(CV_cvodemem,
                    &CV_iout[0],                          /* LENRW   */
                    &CV_iout[1]);                         /* LENIW   */
  CVodeGetIntegratorStats(CV_cvodemem, 
                          &CV_iout[2],                    /* NST     */
                          &CV_iout[3],                    /* NFE     */ 
                          &CV_iout[7],                    /* NSETUPS */ 
                          &CV_iout[4],                    /* NETF    */ 
                          (int *) &CV_iout[8],            /* QU      */
                          (int *) &CV_iout[9],            /* QCUR    */
                          &CV_rout[0],                    /* H0U     */
                          &CV_rout[1],                    /* HU      */ 
                          &CV_rout[2],                    /* HCUR    */ 
                          &CV_rout[3]);                   /* TCUR    */ 
  CVodeGetTolScaleFactor(CV_cvodemem, 
                         &CV_rout[4]);                    /* TOLSFAC */
  CVodeGetNonlinSolvStats(CV_cvodemem,
                          &CV_iout[6],                    /* NNI     */
                          &CV_iout[5]);                   /* NCFN    */
  CVodeGetNumStabLimOrderReds(CV_cvodemem, &CV_iout[10]); /* NOR     */
  
  /* Root finding is on */
  if (CV_nrtfn != 0)
    CVodeGetNumGEvals(CV_cvodemem, &CV_iout[11]);         /* NGE     */
  
  switch(CV_ls) {
  case CV_LS_DENSE:
    CVDenseGetWorkSpace(CV_cvodemem, &CV_iout[12], &CV_iout[13]);  /* LENRWLS,LENIWLS */
    CVDenseGetLastFlag(CV_cvodemem, (int *) &CV_iout[14]);         /* LSTF */
    CVDenseGetNumRhsEvals(CV_cvodemem, &CV_iout[15]);              /* NFELS */
    CVDenseGetNumJacEvals(CV_cvodemem, &CV_iout[16]);              /* NJE */
    break;
  case CV_LS_BAND:
    CVBandGetWorkSpace(CV_cvodemem, &CV_iout[12], &CV_iout[13]);   /* LENRWLS,LENIWLS */
    CVBandGetLastFlag(CV_cvodemem, (int *) &CV_iout[14]);          /* LSTF */
    CVBandGetNumRhsEvals(CV_cvodemem, &CV_iout[15]);               /* NFELS */
    CVBandGetNumJacEvals(CV_cvodemem, &CV_iout[16]);               /* NJE */
    break;
  case CV_LS_DIAG:
    CVDiagGetWorkSpace(CV_cvodemem, &CV_iout[12], &CV_iout[13]);   /* LENRWLS,LENIWLS */
    CVDiagGetLastFlag(CV_cvodemem, (int *) &CV_iout[14]);          /* LSTF */
    CVDiagGetNumRhsEvals(CV_cvodemem, &CV_iout[15]);               /* NFELS */
    break;
  case CV_LS_SPGMR:
  case CV_LS_SPBCG:
  case CV_LS_SPTFQMR:
    CVSpilsGetWorkSpace(CV_cvodemem, &CV_iout[12], &CV_iout[13]);  /* LENRWLS,LENIWLS */
    CVSpilsGetLastFlag(CV_cvodemem, (int *) &CV_iout[14]);         /* LSTF */
    CVSpilsGetNumRhsEvals(CV_cvodemem, &CV_iout[15]);              /* NFELS */
    CVSpilsGetNumJtimesEvals(CV_cvodemem, &CV_iout[16]);           /* NJTV */
    CVSpilsGetNumPrecEvals(CV_cvodemem, &CV_iout[17]);             /* NPE */
    CVSpilsGetNumPrecSolves(CV_cvodemem, &CV_iout[18]);            /* NPS */
    CVSpilsGetNumLinIters(CV_cvodemem, &CV_iout[19]);              /* NLI */
    CVSpilsGetNumConvFails(CV_cvodemem, &CV_iout[20]);             /* NCFL */
  }
}

/***************************************************************************/

void FCV_DKY (realtype *t, int *k, realtype *dky, int *ier)
{
  /* 
     t             is the t value where output is desired
     k             is the derivative order
     F2C_CVODE_vec is the N_Vector containing the solution derivative on return 
  */

  N_VSetArrayPointer(dky, F2C_CVODE_vec);

  *ier = 0;
  *ier = CVodeGetDky(CV_cvodemem, *t, *k, F2C_CVODE_vec);

  N_VSetArrayPointer(NULL, F2C_CVODE_vec);

}

/*************************************************/

void FCV_GETERRWEIGHTS(realtype *eweight, int *ier)
{
  /* Attach user data to vector */
  N_VSetArrayPointer(eweight, F2C_CVODE_vec);

  *ier = 0;
  *ier = CVodeGetErrWeights(CV_cvodemem, F2C_CVODE_vec);

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_CVODE_vec);

  return;
}

/*************************************************/

void FCV_GETESTLOCALERR(realtype *ele, int *ier)
{
  /* Attach user data to vector */
  N_VSetArrayPointer(ele, F2C_CVODE_vec);

  *ier = 0;
  *ier = CVodeGetEstLocalErrors(CV_cvodemem, F2C_CVODE_vec);

  /* Reset data pointers */
  N_VSetArrayPointer(NULL, F2C_CVODE_vec);

  return;
}

/***************************************************************************/

void FCV_FREE ()
{
  CVodeMem cv_mem;

  cv_mem = (CVodeMem) CV_cvodemem;

  free(cv_mem->cv_f_data); cv_mem->cv_f_data = NULL;

  CVodeFree(&CV_cvodemem);

  N_VSetArrayPointer(NULL, F2C_CVODE_vec);
  N_VDestroy(F2C_CVODE_vec);
}

/***************************************************************************/

/* 
 * C function CVf to interface between CVODE and a Fortran subroutine FCVFUN.
 * Addresses of t, y, and ydot are passed to CVFUN, using the
 * routine N_VGetArrayPointer from the NVECTOR module.
 * Auxiliary data is assumed to be communicated by Common. 
 */

int FCVf(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  int ier;
  realtype *ydata, *dydata;
  FCVUserData CV_userdata;

  ydata  = N_VGetArrayPointer(y);
  dydata = N_VGetArrayPointer(ydot);

  CV_userdata = (FCVUserData) f_data;

  FCV_FUN(&t, ydata, dydata, CV_userdata->ipar, CV_userdata->rpar, &ier);

  return(ier);
}
