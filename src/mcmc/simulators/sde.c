#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>

#include "../../headers.h"
#include "../mcmc_headers.h"

/*Globals because I can't be bothered anymore!*/
int postive_definite=1;

void drift2(gsl_vector *mean_cand, gsl_vector *param1, double ddeltat)
{
  double m1,m2,k1,k2,k3,k4,k5,r1,r2;
  k1=VGET(param1,0);
  k2=VGET(param1,1);
  k3=VGET(param1,2);
  k4=VGET(param1,3);
  k5=VGET(param1,4);
  r1=VGET(mean_cand,0);
  r2=VGET(mean_cand,1);
    
  m1=ddeltat*(k1-k3*r1-k5*r1*r2);
  m2=ddeltat*(k2-k4*r2+20*k5*r1*r2);/*Change to 20*k5*r1*r2*/
    
  VSET(mean_cand, 0, m1+VGET(mean_cand, 0));
  VSET(mean_cand, 1, m2+VGET(mean_cand, 1));
}

void diffusion2(gsl_vector *mean_cand, gsl_vector *param1, double ddeltat, gsl_matrix *disp_mat)
{
  double v1,v2,v12,k1,k2,k3,k4,k5,r1,r2;
  k1=VGET(param1, 0);
  k2=VGET(param1, 1);
  k3=VGET(param1, 2);
  k4=VGET(param1, 3);
  k5=VGET(param1, 4);
  r1=VGET(mean_cand, 0);
  r2=VGET(mean_cand, 1);

  v1 = ddeltat*(k1 + k3*r1 + k5*r1*r2);
  v2 = ddeltat*(k2 + k4*r2 + 400*k5*r1*r2);/*Change: 20^2*k5...*/
  v12 = -ddeltat*(20*k5*r1*r2);/*20*k5*...*/
  if((k1+k2+k3+k4) < 0.00000001) {
      postive_definite = 0;
      /*printf("#######################\n");*/
    }
  else
    postive_definite = 1;
  MSET(disp_mat,0,0,v1);
  MSET(disp_mat,1,1,v2);
  MSET(disp_mat,0,1,v12);
  MSET(disp_mat,1,0,v12);
}

void mvn_sample(gsl_vector *mean_cand, gsl_matrix *var)
{
  /* Takes a mean vec, mean and var matrix, 
   * var and gives vector of MVN(mean,var) realisations, x 
   */
  int i, j;
  int dimen = var -> size1;
  double value;
  gsl_matrix *disp;
  gsl_vector *ran;
  gsl_matrix *fast_species;
  
  fast_species = gsl_matrix_alloc(2, 2);
  gsl_matrix_set_identity(fast_species);
  
  for(i=0;i<dimen; i++) {
    if(MGET(var, i, i) <0.00000000001) {
      MSET(var, i, i, 1.0);
      MSET(fast_species, i, i, 0.0);
    }
  }
  
  disp = gsl_matrix_alloc(2, 2);
  ran = gsl_vector_alloc(2);
  gsl_matrix_memcpy(disp, var);
  if(postive_definite == 1) {
    gsl_linalg_cholesky_decomp(disp);
    for(i=0;i<dimen;i++) {
      for (j=i+1;j<dimen;j++) {
        MSET(disp,i,j,0.0);
      }
    }
  }else{
    value = pow(MGET(disp, 0 ,0), 0.5);
    gsl_matrix_set_identity(disp);
    MSET(disp, 0,0, value);
    MSET(disp, 1,1, value);       
  }

  for (j=0;j<dimen;j++) {
    VSET(ran,j,gsl_ran_gaussian(r,1.0));
  }

  /*remove update from slow species*/
  gsl_matrix_mul_elements(disp, fast_species);
    
  /*Add noise to mean cand*/
  gsl_blas_dgemv(CblasNoTrans,1.0, disp, ran, 1.0, mean_cand);
  for(i=0; i<2; i++)  {
    if(VGET(mean_cand,i)<=0.0001 && MGET(fast_species, i, i) > 0.000001)
      VSET(mean_cand,i,0.0001);
  }
  gsl_vector_free(ran);
  gsl_matrix_free(disp);
  gsl_matrix_free(fast_species);
}

double get_min(double t1, double t2) {
  if(t1 < t2) {
    return(t1);
  } else {
    return(t2);
  }
}


/*
  Public functions
*/

/*SDE solution*/
void sde(gsl_vector *par, double tlen, gsl_vector *sps)
{

  double ddeltat, t=0;
  double total_par = 0.0;
  int i;
  for(i=0; i< par->size; i++)
    total_par += VGET(par, i);

  if(total_par < 0.000000000001)
    return;
  gsl_matrix *disp_mat = gsl_matrix_alloc(2, 2);

  while((tlen-t) > 0.000000001) {
    ddeltat = get_min(tlen-t, 0.005);/* SDE time step */
    diffusion2(sps, par, ddeltat, disp_mat);  
    drift2(sps, par, ddeltat);
    mvn_sample(sps, disp_mat);
    t += ddeltat;
  }
  gsl_matrix_free(disp_mat);
}

void SDE(st_part_at* prop_part, double tstep)
{
  if(tstep < 0.0001) return;
  sde(prop_part->params, tstep, prop_part->sps);
}


