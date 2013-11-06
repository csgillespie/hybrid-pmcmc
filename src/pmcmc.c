#include <getopt.h>   

#include "headers.h"
#include "mcmc/mcmc_headers.h"
#include "mcmc/algorithm.h"
#include "mcmc/init.h"
#include "mcmc/generic/commandline.h"
#include "mcmc/simulators/simulators.h"
#include "mcmc/generic/paths.h"
 
/*User Globals
  experiment_file
  * col 1 time, other columns species, not necessary all species
  * Currently assume that if we observe a sp at time t_i, we observe for all time
  prior_file
  * particles determined from length
  * number of columns == number of species. 
  * Order of columns must be the same as experiment_file
  * If prior contains more columns these are assume to be unobserved species
  */

/*Two arguments burnin and thin*/
int main(int argc, char *argv[])
{

  r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (r, 3);   

  int burn = getBurnin(argc, argv);
  int thin = getThin(argc, argv);
  long iters = getIters(argc, argv, burn, thin);
  char *output_file = getOutputFile(argc, argv);
  char *dir = getDir(argc, argv);
  char *simulator_name = getSimulatorName(argc, argv);
  
  dir = createPath(dir, simulator_name);
  dir = createPath(dir, "/");

  st_mcmc_settings *mcmc_settings = initMCMCSettings(burn, thin, iters, dir, output_file); 

  st_model_at *model_at = initModelAttributes();    
  void (*forwardSimulate)(st_part_at* prop_part, double tstep);
  forwardSimulate = getForwardSimulator(argc, argv);

  part(model_at, mcmc_settings, forwardSimulate, dir);  

  return(GSL_SUCCESS);
}
 
