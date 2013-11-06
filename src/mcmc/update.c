#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>


#include "../headers.h"
#include "mcmc_headers.h"
/*Private  functions */

/*particle monte-carlo mean & variance for particle filter*/
double samp_mean(gsl_vector *pars)
{
  int i;
  double smean;
  smean=0.0;

  for (i=0; i<pars->size; i++){
    smean += log(VGET(pars, i));
  }
  smean = smean/(pars->size);
  return(smean);
}


double samp_variance(gsl_vector *pars)
{
  int i;
  double smean, svar;
  svar=0.0;
  smean = samp_mean(pars);

  for(i=0; i<pars->size; i++) {
    svar += (log(VGET(pars, i))-smean)*(log(VGET(pars, i))-smean);
  }
  svar = svar/(pars->size - 1.0);  
  return(svar);
}

void updateParameters(st_part_at *accepted_part, 
                      st_part_at *prop_part,
                      st_mcmc_update *mcmc_update)
{
  int j;
  gsl_vector *prop_pars = prop_part->params;
  gsl_vector *accepted_pars = accepted_part->params;

  gsl_matrix *tuning = mcmc_update->tuning;
  gsl_vector *fixed = mcmc_update->fixed;
  gsl_vector *z = mcmc_update->z;
  double upper[6],lower[6];
  
  lower[0]=-15; lower[1]=-15; lower[2]=-15; lower[3]=-15; lower[4]=-15;  
  upper[0]=15; upper[1]=15; upper[2]=15; upper[3]=15; upper[4]=15;  

  for (j=0; j<prop_pars->size; j++) {
    VSET(z, j, gsl_ran_gaussian(r,1.0));
    VSET(prop_pars, j, log(VGET(accepted_pars, j)));//Log parameters
  }
  gsl_blas_dgemv(CblasNoTrans, 1.0, tuning, z, 1.0, prop_pars);

  prop_part->like = 0.0;
  for(j=0; j<prop_pars->size; j++) {
    if(VGET(fixed, j) > 0.5) {
      VSET(prop_pars, j, VGET(accepted_pars, j));
    } else {
      prop_part->like += log(gsl_ran_flat_pdf(VGET(prop_pars, j),lower[j], upper[j]));
      VSET(prop_pars, j, exp(VGET(prop_pars, j)));
      //      printf("j: %d : %f", j,  VGET(prop_pars, j));
    }
  }  
  //  printf("\n");
}

void updateSpecies(int index, st_parts_at *prior_parts, st_part_at *prop_part)
{           
  gsl_matrix_get_row(prop_part->sps, prior_parts->sps, index);
}

void updateResiduals(int index, st_parts_at *prior_parts, st_part_at *prop_part)
{   
  gsl_matrix_get_row(prop_part->res, prior_parts->res, index);
}

void acceptParameters(st_model_at *model_at, st_part_at *accepted_part,
                      st_part_at *prop_part, int accept)
{
  double u, aprob;
  aprob = prop_part->like - accepted_part->like;
  u = gsl_ran_flat(r, 0.0, 1.0);
  if(log(u) < aprob || accept){
    gsl_vector_memcpy(accepted_part->params, prop_part->params);
    accepted_part->like = prop_part->like;
  }
}
             

void switchParticles(st_parts_at *prior_parts, st_parts_at *accepted_parts)
{
  gsl_matrix_memcpy(prior_parts->sps, accepted_parts->sps);
  gsl_matrix_memcpy(prior_parts->res, accepted_parts->res);    
}  

void updateParticles(int iters, st_part_at *accepted_part, st_parts_at *accepted_parts)
{                
  gsl_matrix_set_row(accepted_parts->sps, iters, accepted_part->sps);
  gsl_matrix_set_row(accepted_parts->res, iters, accepted_part->res);            
    
}       




