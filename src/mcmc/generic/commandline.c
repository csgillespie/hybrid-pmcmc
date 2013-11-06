#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <getopt.h>  

#include "../mcmc_headers.h"
#include "../simulators/simulators.h"
#include "commandline.h"

char *getSimulatorOptions()
{
  char *options = "b:s:t:d:n:";
  return(options);
}

char *getTimingOptions()
{
  char *options = "b:s:t:d:";
  return(options);
}


/* Public functions */
long getIters(int argc, char *argv[], int burn, int thin) 
{
  int opt;
  long iters=0;
  char *options = getSimulatorOptions();
  while ((opt = getopt (argc, argv, options)) != -1) {
    switch(opt)
    {
    case 'n' :
      iters = atoi(optarg);
      break;
    }
  }
  if(iters == 0) {
    iters = burn + thin;
  } else { 
    iters = burn + thin*iters;
  }
  optind = 1;
  return(iters);
}

char *getDir(int argc, char *argv[]) 
{
  int opt;
  char *dir = NULL;
  char *options = getSimulatorOptions();
  while ((opt = getopt (argc, argv, options)) != -1) {
    switch(opt)
    {
    case 'd' :
      dir = optarg;
      break;
    }
  }
   if(dir == NULL) {
    printf("Directory not set. Using current dir\n");
    dir = "";
    //    exit(GSL_FAILURE);
  }
  optind = 1;
  return(dir);
}




char *getSimulatorName(int argc, char *argv[]) 
{
  int opt;
  char *output_file = NULL;
  char *options = getSimulatorOptions();

  while ((opt = getopt (argc, argv, options)) != -1) {
    switch(opt)
    {
    case 's' :
      if(strcmp(optarg, "gillespie")==0){
        output_file = "gillespie";
      } else if (strcmp(optarg, "hybridSDE")==0){
        output_file = "hybridSDE";
      } else if (strcmp(optarg, "hybridLNA")==0){
        output_file = "hybridLNA";
      } else {
        printf("Error: Simulator not recognised\n");
        exit(GSL_FAILURE);
      }
      break;
    }
  }
  if(output_file == NULL) {
    printf("Error: Simulator not set\n");
    exit(GSL_FAILURE);
  }
  optind = 1;
  return(output_file);
}


char *getOutputFile(int argc, char *argv[]) 
{
  int opt;
  char *output_file = NULL;
  char *options = getSimulatorOptions();

  while ((opt = getopt (argc, argv, options)) != -1) {
    switch(opt)
    {
    case 's' :
      if(strcmp(optarg, "gillespie")==0){
        output_file = "gillespie.csv";
      } else if (strcmp(optarg, "hybridSDE")==0){
        output_file = "hybridSDE.csv";
      } else if (strcmp(optarg, "hybridLNA")==0){
        output_file = "hybridLNA.csv";
      } else {
        printf("Error: Simulator not recognised\n");
        exit(GSL_FAILURE);
      }
      break;
    }
  }
  if(output_file == NULL) {
    printf("Error: Simulator not set\n");
    exit(GSL_FAILURE);
  }
  optind = 1;
  return(output_file);
}

void *getForwardSimulator(int argc, char *argv[]) 
{
  int opt;
  void (*forwardSimulate)(st_part_at* prop_part, double tstep) = NULL;
  char *options = getSimulatorOptions();
  while ((opt = getopt (argc, argv, options)) != -1) {
    switch(opt)
      {
      case 's' :
        if(strcmp(optarg, "gillespie")==0){
          forwardSimulate = gillespie;
        } else if (strcmp(optarg, "hybridSDE")==0){
          forwardSimulate = hybridSDE;
        } else if (strcmp(optarg, "hybridLNA")==0){
          forwardSimulate = hybridLNA;
        } else {
          printf("Error: need to set a simulator\n");
          exit(GSL_FAILURE);
        }
        break;
      }
  }
  if(forwardSimulate == NULL) {
    printf("Error: Simulator not set\n");
    exit(GSL_FAILURE);
  }
  optind = 1;
  return(forwardSimulate);
}




int getBurnin(int argc, char *argv[]) 
{

  int opt, burnin=0;
  char *options = getSimulatorOptions();

  while ((opt = getopt (argc, argv, options)) != -1) {
    switch(opt) 
      {
      case 'b':
        burnin = atoi(optarg);
        break;
      }
  }
  optind = 1;
  return(burnin);

}

int getThin(int argc, char *argv[]) 
{
  int opt, thin=1;
  char *options = getSimulatorOptions();

  while ((opt = getopt (argc, argv, options)) != -1) { 
    switch(opt) 
      {
      case 't':
        thin = atoi(optarg);
        break;
      }
  }
  optind = 1;
  return(thin);

}
