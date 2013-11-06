#include <math.h>
#include <gsl/gsl_randist.h>

#include "../../headers.h"
#include "../mcmc_headers.h"

void gillespie(st_part_at* prop_part, double tstep)
{
    
  double alpha_x, mu_x, beta, alpha_y, mu_y, u, time;
  double X, Y, rate;
  X = VGET(prop_part->sps, 0);
  Y = VGET(prop_part->sps, 1);

  
  alpha_x=VGET(prop_part->params, 0);
  alpha_y=VGET(prop_part->params, 1);
  mu_x=VGET(prop_part->params, 2);
  mu_y=VGET(prop_part->params, 3);
  beta=VGET(prop_part->params, 4);

  time = 0;

  while(time < tstep) {

    rate = alpha_x + X*mu_x + beta*X*Y + Y*mu_y + alpha_y;
        
    u = gsl_ran_flat(r, 0.0, 1.0);
    time += -log(u)/rate;
    if(time > tstep)
      break;
    u = gsl_ran_flat(r, 0.0, 1.0);

    if(u<alpha_x/rate) {
      X += 1;
    } else if(u< (alpha_x+X*mu_x)/rate) {
      X -= 1;
    } else if(u <(alpha_x+X*mu_x + beta*X*Y)/rate) {
      X -= 1;
      Y += 20;
    } else if(u <(alpha_x + X*mu_x + beta*X*Y + Y*mu_y)/rate) {
      Y -= 1;
    } else {
      Y += 1;
    }
  }
   
  VSET(prop_part->sps, 0, X);
  VSET(prop_part->sps, 1, Y);
}


