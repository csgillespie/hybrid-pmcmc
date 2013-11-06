#include "../../headers.h"
#include "../mcmc_headers.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <time.h>
#include <sys/timeb.h>


#include "Generic_c/LinearAlgebra/cslinalg.h"
#include "Generic_c/Print/arrvec.h"
#include "LNA_generic_hybrid/struct_generic.h"
#include "LNA_generic_hybrid/read.h"
#include "LNA_generic_hybrid/genericODE.h"
#include "LNA_generic_hybrid/system_spec.h"
#include "LNA_generic_hybrid/random.h"
#include "hybrid_lna.h"

#define VSET gsl_vector_set
#define VGET gsl_vector_get
#define MSET gsl_matrix_set
#define MGET gsl_matrix_get


//#define PRINT_STATEMENTS


// ***************************

  



void hybrid_lna(double t_interval, double *initial, double *pars) {
  SDEInfo *pSI = (SDEInfo*)pSIv;
  ReactionInfo *pRI=pSI->pRI;

  int i = 111111;
  double this_int1,this_int2;
  double *pzero; // S*S*R array of zeros
  int S=pSI->S;
  int R=pSI->R;
  //  gsl_rng_set(r,seed);

  pzero=(double*)malloc(S*S*R*sizeof(double));
  for (i=0;i<S*S*R;i++) {
    pzero[i]=0.0;
  }
  memcpy(pRI->pc,pars,R*sizeof(double));

  ODEIntegrateSimulateOne(t_interval, pSI, pODEIv, pzero, initial, r,
			  &this_int1,&this_int2);
 

  free(pzero);

} 

void hybridLNA(st_part_at* prop_part, double tstep) {

  double *sps = prop_part->sps->data;
  double *pars = prop_part->params->data;
  hybrid_lna(tstep, sps, pars);
}


void *SetUpSIRI(void) {
  ReactionInfo *pRI=(ReactionInfo*)malloc(sizeof(ReactionInfo));
  SDEInfo *pSI=(SDEInfo*)malloc(sizeof(SDEInfo));

  AllocReactionInfo(pRI);
  AllocSDEInfo(pRI,pSI);
  return (void*)pSI;
}

void DestroySIRI(void *pSIv) {
  SDEInfo *pSI=(SDEInfo*)pSIv;
  ReactionInfo *pRI=pSI->pRI;

  DeAllocSDEInfo(pSI);
  free(pSI);
  DeAllocReactionInfo(pRI);
  free(pRI);
}

void *pgetvoidRI(void *pSIv) {
  SDEInfo *pSI=(SDEInfo*)pSIv;
  return (void*) pSI->pRI;
}
 

double test() {return(1.0);}
/* int main(int argc, char *argv[]) {
 *   void *pODEIv, *pSIv;
 *   double t_interval;
 *   double sc=5;
 *   double pars[5];
 *   double initial[2];
 *   
 *   // Set up
 *   
 *   pSIv=SetUpSIRI();
 *   pODEIv = pODEInfoAlloc(pgetvoidRI(pSIv));
 * 
 * 
 *   // COLIN: change from here ...
 *   t_interval=1000;
 *   sc=5;
 *   // GG - these are the values in your email of 18/11/2012
 *   pars[0]=2; pars[1]=sc; pars[2]=0.02; pars[3]=1; pars[4]=0.02/sc;
 *   initial[0]=0; initial[1]=0;
 * 
 * 
 *   
 *   hybridLNA(t_interval, initial, pars, pSIv, pODEIv);
 * 
 *   // print the simulated values 
 *   printdvec(2,initial);
 * 
 *   // ... to here
 * 
 *   // Free
 * 
 *   ODEInfoDeAlloc(pODEIv);
 *   DestroySIRI(pSIv);
 * 
 *   return(0);
 * } */

