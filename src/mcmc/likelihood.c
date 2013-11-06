
#include "../headers.h"
#include "mcmc_headers.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
double calculatePoissonPerSps(double observed, double predicted) {

    double ll;    
    observed = round(observed); predicted = round(predicted);

    if(predicted > 0.5) {
      ll = log(gsl_ran_poisson_pdf(observed, predicted)); 
    } else {
      ll = observed*log(0.1) + (1.0-predicted)*log(0.9);
      if(observed > 1.0) {
        ll = GSL_NEGINF;
      }     
    }

    return(ll);
}  

/*Public functions*/
double calculatePoissonLike(double *observed, st_part_at *prop_part) {

    gsl_vector *predicted = prop_part->sps;
    double ll, pred, obs;
    int i;
    ll = 0.0;

    for(i=0; i<predicted->size; i++) {
      obs = observed[i];
      pred = VGET(predicted, i);
      ll += calculatePoissonPerSps(obs, pred);
    }

    return(exp(ll));
}














