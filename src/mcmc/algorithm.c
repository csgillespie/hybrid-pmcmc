#include <math.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <sys/timeb.h>


#include "../headers.h"
#include "mcmc_headers.h"
#include "init.h"
#include "likelihood.h"
#include "Output.h"
#include "update.h"
#include "simulators/simulators.h"

#include "generic/timings.h"


void part(st_model_at *model_at, 
          st_mcmc_settings *mcmc_settings,
          void (*forwardSimulate)(st_part_at* prop_part, double tstep),
          char *dir)
{
  
 /* For timings */
  struct timeb t_start, t_end;                          


  int k, i, index=0, T;
  double *wts;
  double partial_like;
  int thin = mcmc_settings->thin;
  int burn = mcmc_settings->burn;

  gsl_ran_discrete_t *g;
  st_part_at *accepted_part, *prop_part;
  st_parts_at *accepted_parts, *prior_parts, *unsampled_parts;
  st_mcmc_update *mcmc_update;

  st_data *data_at = initData(model_at, dir);
  st_data *root_data = data_at;

  prior_parts = initPrior(dir);
  unsampled_parts = initAcceptedParts(prior_parts);
  accepted_parts = initAcceptedParts(prior_parts);

  accepted_part = initPart(prior_parts, dir);
  prop_part = initPart(prior_parts, dir);

  mcmc_update = initMCMCUpdate(dir);
  model_at->no_params = accepted_part->params->size;
  model_at->no_sps = prior_parts->sps->size2;

  mcmc_settings->parts = prior_parts->sps->size1;
  outputMCMCSettings(mcmc_settings);
  outputPartSettings(accepted_part);
  outputPartsSettings(accepted_parts);

  /*Used for the hybrid LNA*/
  pSIv = SetUpSIRI();
  pODEIv = pODEInfoAlloc(pgetvoidRI(pSIv));

  ftime(&t_start);
  wts = calloc(sizeof(double), mcmc_settings->parts);
  for(i=0; i<mcmc_settings->iters; i++){
    updateParameters(accepted_part, prop_part, mcmc_update);
    T = 0;
    data_at = root_data;    
    switchParticles(accepted_parts, prior_parts);
        
    while(data_at != NULL) {
      //outputData(data_at, model_at->no_obs_sps);
      partial_like = 0.0; T++;
      for(k=0; k<(mcmc_settings->parts); k++) {
        updateSpecies(k, accepted_parts, prop_part); 
        updateResiduals(k, accepted_parts, prop_part);
        /*forward simulate*/
        forwardSimulate(prop_part, data_at->tstep);
        /* Calculate likelihood */
        wts[k] = calculatePoissonLike(data_at->obser, prop_part);
                
        partial_like += wts[k]; 
        updateParticles(k, prop_part, unsampled_parts);
      }
      prop_part->like += log(partial_like);
      
      g = gsl_ran_discrete_preproc(mcmc_settings->parts, (const double *) wts);
      for(k=0; k<mcmc_settings->parts; k++){
        index = gsl_ran_discrete(r, g);
        gsl_matrix_get_row(prop_part->sps, unsampled_parts->sps, index);
        gsl_matrix_get_row(prop_part->res, unsampled_parts->res, index);
        updateParticles(k, prop_part, accepted_parts);
      }
      gsl_ran_discrete_free(g);

      data_at = data_at->next;
    } 

    prop_part->like -= T*log(mcmc_settings->parts);
    // outputParameterVector(accepted_part);
    //outputParameterVector(prop_part);
    acceptParameters(model_at, accepted_part, prop_part, mcmc_update->always_accept);
    //  printf("###i=%d\n",i );
    if(((i-burn) >= 0) && ((i-burn) % thin == 0)) {
       outputParameters(i, mcmc_settings, accepted_part);
    }
  }
  ftime(&t_end);
  outputTimings(toMilliseconds(t_start, t_end), mcmc_settings);
  /*Hybrid_LNA*/
  ODEInfoDeAlloc(pODEIv);
  DestroySIRI(pSIv);
  gsl_rng_free(r);
  free(wts);
  //return(GSL_SUCCESS);
}



