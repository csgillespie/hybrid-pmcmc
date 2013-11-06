#include <gsl/gsl_randist.h>
#include <math.h>

#include "../../headers.h"
#include "../mcmc_headers.h"
#include "sde.h"


/*****************
Private functions
*******************/
double min(double value1, double value2)
{
  if(value1 > value2)
    return(value2);
  return(value1);
}

double max(double value1, double value2)
{
  if(value1 < value2)
    return(value2);
  return(value1);
}

double round(double value1)
{
  if(value1 > 0.5)
    return(max(1, value1));
  return(0.0);
}

/*
  Hybrid method
*/
void updateHazard(gsl_vector *haz, gsl_vector *params, gsl_vector *sps)
{
  VSET(haz, 0, VGET(params, 0));
  VSET(haz, 1, VGET(params, 1));
  VSET(haz, 2, VGET(params, 2)*VGET(sps, 0));        
  VSET(haz, 3, VGET(params, 3)*VGET(sps, 1));        
  VSET(haz, 4, VGET(params, 4)*VGET(sps, 0)*VGET(sps, 1)); 
}           

/*Corresponds to equation 10*/
int checkSpecies(gsl_matrix *PostPre, gsl_vector *sps, int reaction_no)
{
  int i;
  double cond, epsilon_st, N_st;

  N_st = 15; epsilon_st = 0.25;
  for(i=0; i<sps->size; i++) {
    cond = fabs(MGET(PostPre, reaction_no, i))*N_st;
    if(cond > (epsilon_st*VGET(sps, i))) return(0);
  }
  return(1);
}


/*A bit inefficient to have two separate functions. Just easier
 * Corresponds to equation 9
 */
int checkHazard(gsl_vector *haz,gsl_matrix *PostPre, gsl_vector *sps, 
                int reaction_no, double deltat) 
{
  int i;
  double cond,  epsilon=0.25;

  for(i=0; i<sps->size; i++) {
    cond =  fabs(MGET(PostPre, reaction_no, i));
    cond *= max(1.0, VGET(haz, reaction_no)*deltat);
    if(cond > epsilon*VGET(sps, i)) return(0);
  }
  return(1);
}

int running_totals1[5] = {0, 0, 0, 0, 0};
double checkReactions(gsl_vector *haz, gsl_vector *residuals, gsl_vector *fast_params, 
                      double deltat,  gsl_vector *sps, gsl_matrix *PostPre)
{
  int i;
  double small_residual = -1000000.0;

  double total_hazard = 0.0;
  for(i=0; i < fast_params->size; i++) {
    if(checkHazard(haz, PostPre, sps, i, deltat) && 
       checkSpecies(PostPre, sps, i))  {
      VSET(residuals, i, small_residual);
    } else if(VGET(residuals, i) < (small_residual + 10000))  {
      total_hazard += VGET(haz, i);
      VSET(residuals, i, log(gsl_ran_flat(r, 0.0, 1.0)));
      VSET(fast_params, i, 0.0);
      //      running_totals1[i]++;
    } else {
      total_hazard += VGET(haz, i);
      VSET(fast_params, i, 0.0);
      //      running_totals1[i]++;
    }
    //printf("%d ", running_totals1[i]);
  }
  //  printf("\n");
  return total_hazard;        
}

void hybridSim(st_part_at* prop_part,
               double maxtime, 
               void (*forwardsimulate)(gsl_vector *params, double maxtime, gsl_vector* sps)/*a pointer to the forward simulator method*/
               ) {
    
  double t,  tau=0;
  int zc=-1, i, mu=-1;
  double  deltat, Deltat=0.1;
  deltat = Deltat;
    
  gsl_vector *sps = prop_part->sps;
  gsl_vector *params = prop_part->params;
  gsl_vector *residuals = prop_part->res;        
    
  gsl_vector *fast_params;
  fast_params = gsl_vector_alloc(params->size);
  gsl_vector_memcpy(fast_params, params);
    
  gsl_matrix *PostPre;
  PostPre = gsl_matrix_alloc(5, 2);
  gsl_matrix_set_zero(PostPre);
  MSET(PostPre, 0, 0, 1);
  MSET(PostPre, 1, 1, 1);    
  MSET(PostPre, 2, 0, -1);        
  MSET(PostPre, 3, 1, -1);            
  MSET(PostPre, 4, 0, -1);                
  MSET(PostPre, 4, 1, 20);/*Change to 1-> 20*/               
    
  gsl_vector *haz;
  haz = gsl_vector_alloc(5);
  gsl_vector_set_zero(haz);
  double toutput = 0.0;
  t = 0.0;
  while (t<maxtime){
    deltat = min(deltat, maxtime - t );
    gsl_vector_memcpy(fast_params, params);
    //    printf("sps0=%f, sps1=%f\n", VGET(sps, 0),VGET(sps, 1));
    /*Update Hazard function*/
    updateHazard(haz, params, sps);
    checkReactions(haz, residuals, fast_params, deltat, sps, PostPre);
        
    /*Evaluate residuals of jump equations*/
    gsl_vector_scale(haz, deltat);
    gsl_vector_add(residuals, haz);
    zc = 0;
    for(i=0; i<5; i++) {
      if(VGET(residuals, i)>=0.0) {
        zc++;
        mu = i;
      }
    }
    if(zc == 0) {
      forwardsimulate(fast_params, deltat, sps);
      t += deltat;
      deltat = Deltat;
    }
        
    if(zc == 1) {
      /*Generate time to reaction*/
      tau = -(VGET(residuals, mu) - VGET(haz, mu))/(VGET(haz, mu)/deltat);
      gsl_vector_sub(residuals, haz);
      gsl_vector_scale(haz, tau/deltat);
      gsl_vector_add(residuals, haz);
      VSET(residuals, mu, log(gsl_ran_flat(r, 0.0, 1.0)));
      deltat = tau;
      forwardsimulate(fast_params, deltat, sps);
            
      /*Add on effect of reaction mu*/
      VSET(sps, 0, round(VGET(sps, 0) + MGET(PostPre, mu, 0)));
      VSET(sps, 1, round(VGET(sps, 1) + MGET(PostPre, mu, 1)));            
      t += deltat;
      deltat = Deltat;
    }
    
    if(zc > 1) {
      gsl_vector_sub(residuals, haz);
      deltat = deltat/5.0;
    }

    if(zc < 2 && t > toutput) {
      /* printf("%f,%f,%f,%f,%f,%f\n", 
       *        toutput, VGET(fast_params, 0),
       *        VGET(fast_params, 1), VGET(fast_params, 2),
       *        VGET(fast_params, 3),VGET(fast_params, 4)); */
        toutput += 1.0;
        }
  }

  gsl_vector_free(haz);
  gsl_vector_free(fast_params);
  gsl_matrix_free(PostPre);
}


/*Public Methods*/

/* void hybridODE(st_part_at* prop_part, double maxtime)
 * {
 *   hybridSim(prop_part, maxtime, ode);
 * }     */


void hybridSDE(st_part_at* prop_part, double maxtime)
{
  hybridSim(prop_part, maxtime, sde);
}    





















