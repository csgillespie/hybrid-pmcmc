#include "headers.h"
#include "mcmc/mcmc_headers.h"
#include "mcmc/simulators/simulators.h"
#include "mcmc/generic/timings.h"
#include "mcmc/generic/paths.h"
#include "mcmc/generic/commandline.h"
void *pODEInfoAlloc(void *pRIv);
void ODEInfoDeAlloc(void *pOI);

#include <time.h>
#include <sys/timeb.h>
#include <getopt.h>
#include <gsl/gsl_randist.h>
#include <math.h>

st_part_at * timingsInitPrior(double X, double Y)
{
  st_part_at *part;
  part = (st_part_at *) malloc(sizeof(st_part_at));
  part->params = gsl_vector_alloc(5);
  part->sps = gsl_vector_alloc(2);
  part->res = gsl_vector_alloc(5);

  VSET(part->sps, 0, X);
  VSET(part->sps, 1, Y);
  return(part);
}    

void updateRes(gsl_vector *res) 
{
  int i;
  for(i=0; i<res->size; i++)  
    VSET(res, i, log( gsl_ran_flat(r, 0.0, 1.0)));
}

int main(int argc, char *argv[])
{
  /* Command line arguments
   * -s simulator
   * -n no of simulations
   */
  void (*forwardSimulate)(st_part_at* prop_part, double tstep);
  forwardSimulate = getForwardSimulator(argc, argv);
  int n = getIters(argc, argv, 0, 1);

  char *output_file = getOutputFile(argc, argv);
  char *dir = getDir(argc, argv);
  char *path = createPath(dir, output_file);
  printf("Writing to %s\n", path);
  FILE *out = fopen(path, "w");
  int i, j, k;
  r = gsl_rng_alloc(gsl_rng_mt19937);
  
  /* For timings */

  
  /* Values for scaling parameter*/
  int N_C = 4;

  double sc;
  /*The following three vectors were generated in R and copied in */
  double sc_vec[41] =  {1,10, 100, 1000, 10000};
  double X1[41] = {0,0, 0, 0, 0};
  double X2[41] = {0, 0, 0, 0, 0};
  
  /* Initialise LNA */
  pSIv = SetUpSIRI();
  pODEIv = pODEInfoAlloc(pgetvoidRI(pSIv));
  
  /* Initialise parameters
   * Other paramters and IC are changed in the 
   * loop
   */
  st_part_at *part;
  part = timingsInitPrior(0, 0);
  VSET(part->params, 0, 2.0); 
  VSET(part->params, 2, 1.0/50.0); 
  VSET(part->params, 3, 1.0);
  
  
  for(k=0; k<N_C; k++) {
    sc = sc_vec[k];
    VSET(part->params, 1,sc); VSET(part->params, 4, 1.0/(50*sc));
    for(i=0; i<n; i++){
       VSET(part->sps, 0, X1[k]); VSET(part->sps, 1, X2[k]);
      //VSET(part->sps, 0, 0); VSET(part->sps, 1, 0);
       for(j=0; j < 50; j++){
         fprintf(out, "%d,%d,%d,%d\n", (int) sc, j, 
                 (int) VGET(part->sps, 0),(int) VGET(part->sps, 1));
         updateRes(part->res);
         forwardSimulate(part, 1);  
       }
       fprintf(out, "%d,%d,%d,%d\n", (int) sc, j, 
               (int) VGET(part->sps, 0), (int) VGET(part->sps, 1));
       printf("%d\n", i);
       
    }
  }
  
  ODEInfoDeAlloc(pODEIv);
  DestroySIRI(pSIv);
  gsl_rng_free(r);
  return(GSL_SUCCESS);
}
