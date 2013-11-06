#include "VG_integrate.h"
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h> // added by CS

/* Handle errors within R */
//#include <R.h> - CS: not using R

extern void VG_integrate(double *y,double *t,double * tout,
    double * theta,int  * odesize,
    void (*devis)(), double * reltol, double * abstol)
{
  int lwr= (*odesize) * (*odesize) +22 + (*odesize) * 16;
  int liw =(*odesize)+20;
  double          rwork[lwr];
  int             iwork[liw];
  double  atol,  rtol;
  int itol = 1; //Scalar tolerances
  int  itask=1, istate=1, iopt=0, jt=2, i;

  iwork[5]=8000; // does VG mean iwork[6]?

  rtol= *reltol;
  atol= *abstol;

  lsoda_(devis, odesize, y, t, tout, &itol, &rtol, &atol, 
      &itask, &istate, &iopt, rwork, &lwr, iwork, &liw, NULL, 
      &jt,  theta);

  assert(istate >0);
  //  if (istate <= 0) {
  //  error("error istate = %d\n", istate);
  //}
}
