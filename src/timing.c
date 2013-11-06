#include "headers.h"
#include "mcmc/mcmc_headers.h"
#include "mcmc/simulators/simulators.h"


#include "mcmc/generic/timings.h"
#include "mcmc/generic/paths.h"
#include "mcmc/generic/commandline.h"
#include <time.h>
#include <sys/timeb.h>
#include <getopt.h>
#include <gsl/gsl_randist.h>


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
  int i, k;
  r = gsl_rng_alloc(gsl_rng_mt19937);
  
  /* For timings */
  struct timeb t_start, t_end;                          
  double t_diff;  
  
  /* Values for scaling parameter*/
  int N_C = 41;
  double sim_time;/*Simulation time in secs, for a single sim*/
  double sc;
  /*The following three vectors were generated in R and copied in */
  double sc_vec[41] =  {1,1.3,1.6,2,2.5,3.2,4,5,6.3,7.9,10, 
                        12.6,15.8,20,25.1,31.6,39.8,50.1, 63.1,79.4,100, 
                        125.9,158.5,199.5,251.2,316.2,398.1,501.2,631,794.3,
                        1000,1258.9,1584.9,1995.3,2511.9,3162.3,3981.1,5011.9,6309.6,7943.3,
                        10000};

  double X1[41] = {1,1,1,2,2,3,4,5,6,7,9,11,12,14,16,17,19,20,
                   21,21,22,23,23,23,24,24,24, 24,24,24,24,24,24,24,24,24,24,24,24,24,24};
  
  double X2[41] = {20,20,20,20,20,20,20,20,20,20,21,21,22,
                   24,26,28,32,36,43,50,60,73,89,110,135,168,209,260,325,407,
                   510,639,802,1007,1265,1591,2000,2515,3164,3981,5010};
  
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
    ftime(&t_start);
    for(i=0; i<n; i++){
       VSET(part->sps, 0, X1[k]); VSET(part->sps, 1, X2[k]);
      //VSET(part->sps, 0, 0); VSET(part->sps, 1, 0);
      updateRes(part->res);
      forwardSimulate(part, 100);  
    }
    ftime(&t_end);
    t_diff = toMilliseconds(t_start, t_end);
    sim_time = (t_diff/1000.0)/n;
    printf("%f, %f\n", sc, sim_time);
    fprintf(out, "%f, %f\n", sc, sim_time);
  }
  
  ODEInfoDeAlloc(pODEIv);
  DestroySIRI(pSIv);
  gsl_rng_free(r);
  return(GSL_SUCCESS);
}
