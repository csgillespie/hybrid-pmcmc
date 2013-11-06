#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "../headers.h"
#include "mcmc_headers.h"
#include "generic/paths.h"



gsl_vector *readVector(char *filename)
{
  int line_length = 500;    
  int ncols;
       
  FILE* f;  
  char *pch;
  char line[line_length];
  gsl_vector *particles;
    
  f=fopen(filename, "r");
  if(NULL==f) {
    fprintf(stderr, "Cannot open file %s\n", filename);
    exit(1);
  }

  ncols = 0;
  /*Scan once to get the dimensions
    there doesn't seem to be a realloc matrix function
    Set pch = fgets to avoid error. Really should check
  */
  pch = fgets(line, line_length, f);
  pch = strtok(line,",");
  while(pch != NULL ) {
    ncols++;
    pch = strtok(NULL,",");
  }
  fclose(f);
        
  /*Create vector and fill up*/
  particles = gsl_vector_alloc(ncols);
  ncols = 0;
  f=fopen(filename, "r");
  pch = fgets(line, line_length, f); 
  pch = strtok(line,",");
  while(pch != NULL ) {
    VSET(particles, ncols, atof(pch));
    ncols++;
    pch = strtok(NULL,",");
  }
    
  fclose(f);
  return(particles);    
}    

gsl_matrix *readMatrix(char *filename)
{
  int line_length = 500;    
  int nrows, ncols;
       
  FILE* f;  
  char *pch;
  char line[line_length];
  gsl_matrix *particles;
    
  f = fopen(filename, "r");
  if(NULL==f) {
    fprintf(stderr, "Cannot open file %s\n", filename);
    exit(1);
  }
  nrows = 0; ncols = 0;
  /*Scan once to get the dimensions
    there doesn't seem to be a realloc matrix function
  */

  while(fgets(line, line_length, f) != NULL){
    pch = strtok(line,",");
    while(nrows == 0 && pch != NULL ) {
      ncols++;
      pch = strtok(NULL,",");
    }
    nrows++;
  }
  
  fclose(f);
        
  /*Create matrix and fill up*/
  particles = gsl_matrix_alloc(nrows, ncols);
  nrows = 0; ncols = 0;
  f=fopen(filename, "r");
      
  while(fgets(line, line_length, f) != NULL){
    pch = strtok(line,",");
    while(pch != NULL ) {
      MSET(particles, nrows, ncols, atof(pch));
      ncols++;
      pch = strtok(NULL,",");
    }
    ncols = 0;
    nrows++;
  }
  fclose(f);
    
  return(particles);    
}    

/* Public initialisations functions */
st_parts_at * initPrior(char* dir)
{
  dir = addInputPath(dir);
  st_parts_at *prior_parts;
  char *residual_file= createPath(dir, "residuals.csv");
  char *species_prior_file= createPath(dir, "species_prior.csv");

  prior_parts = (st_parts_at *) malloc(sizeof(st_part_at));
  prior_parts->sps = readMatrix(species_prior_file);
  prior_parts->res = readMatrix(residual_file);
    
  return(prior_parts);
}    


/* MCMC Var-Cov tuning matrix */
st_mcmc_update *initMCMCUpdate(char *dir)
{
  st_mcmc_update *mcmc_update;

  dir = addInputPath(dir);
  char *tuning_file = createPath(dir, "tuning.csv");
  gsl_matrix *tuning = readMatrix(tuning_file);

  char *fixed_file = createPath(dir, "fixed.csv");
  gsl_vector *fixed = readVector(fixed_file);

  if(tuning->size1 != tuning->size2 || tuning->size1 != fixed->size) {
    printf("Error in initMCMCUpdate\n");
    exit(GSL_FAILURE);
  }
  gsl_vector *z = gsl_vector_alloc(fixed->size);

  mcmc_update = (st_mcmc_update *) malloc(sizeof(st_mcmc_update));
  mcmc_update->tuning = tuning;
  mcmc_update->fixed=fixed;
  mcmc_update->z = z;
  mcmc_update->always_accept = 1;

  int i;
  for(i=0; i<fixed->size; i++) 
    mcmc_update->always_accept *= (int) VGET(mcmc_update->fixed, i);

  return(mcmc_update);
}


st_part_at *initPart(st_parts_at *prior_parts, char *dir)
{
  dir = addInputPath(dir);
  st_part_at *part;
  char *parameter_file = createPath(dir, "parameters_prior.csv");

  part = (st_part_at *) malloc(sizeof(st_part_at));
  part->params = readVector(parameter_file);
  part->sps = gsl_vector_alloc(prior_parts->sps->size2);
  part->res = gsl_vector_alloc(prior_parts->res->size2);
  part->like = -100000000.0;/*Set low log likelihood*/
    
  return(part);
}

st_parts_at *initAcceptedParts(st_parts_at *prior_parts)
{
  int no_parts, no_sps;
    
  no_parts = prior_parts->sps->size1;
  no_sps = prior_parts->sps->size2;
    
  st_parts_at *parts;
  parts = (st_parts_at *) malloc(sizeof(st_part_at));
  parts->sps = gsl_matrix_alloc(no_parts, no_sps);
  parts->res = gsl_matrix_alloc(no_parts, prior_parts->res->size2);

  return(parts);
}


st_mcmc_settings * initMCMCSettings(int burn, int thin, long iters,
                                    char *dir, char *output_file)
{
  st_mcmc_settings *mcmc_settings;
  char *filename;

    
  mcmc_settings = (st_mcmc_settings *) malloc(sizeof(st_mcmc_settings));
  mcmc_settings->output_file = (FILE *) malloc(sizeof(FILE));
  mcmc_settings->timing_file = (FILE *) malloc(sizeof(FILE));
  mcmc_settings->burn = burn;
  mcmc_settings->thin = thin;
  mcmc_settings->iters = iters;


  filename = createPath("../output/", dir);
  filename = createPath(filename, output_file);
  FILE* out = fopen(filename, "w");
  if(NULL==out) {
    printf("Cannot open data file %s\n", filename);
    exit(1);
  }

  mcmc_settings->output_file = out;

  filename = createPath("../output/", dir);
  filename = createPath(filename, "timing.csv");
  FILE* time_out = fopen(filename, "w");
  if(NULL==time_out) {
    printf("Cannot open data file %s\n", filename);
    exit(1);
  }
  mcmc_settings->timing_file = time_out;


  printf("###########Output files###################\n");
  printf("%s\n", filename);
  free(filename);
  return(mcmc_settings);
}

st_model_at *initModelAttributes()
{
  st_model_at * model_at;
  model_at = (st_model_at *) malloc(sizeof(st_model_at));
  return(model_at);
}

    
st_data * initData(st_model_at *model_at, char *dir)
{
  dir = addInputPath(dir);
  int i,j;
  int line_length = 80;
    
  double old_tstep = 0.0;
  FILE* f;  
  char *pch;

  char line[line_length];
  st_data *current_data,  *root_data;
  char *experiment_file= createPath(dir, "full.csv");
   
  root_data = (st_data *) malloc(sizeof(st_data));
  root_data->obser = (double *) malloc(sizeof(double));
  root_data->next = NULL;
  current_data = root_data;
  f = fopen(experiment_file, "r");
  if(NULL == f) {
    fprintf(stderr, "Cannot open data file %s\n", experiment_file);
    exit(1);
  }

  i=0;
  j=0;
  while(fgets(line, line_length, f) != NULL){
    j=0;
    if(i>0){
      current_data->next = (st_data *) malloc(sizeof(st_data ));
      current_data = (current_data->next);
      current_data->obser = (double *) malloc(sizeof(double));
    }
    pch = strtok(line,",");
       
    while (pch != NULL ) {
      if(j == 0) {
        current_data->tstep = atof(pch) - old_tstep;
        old_tstep = atof(pch);
      } else {
        current_data->obser = (double *) realloc(current_data->obser, (j*sizeof(double)));
        current_data->obser[j-1] = atof(pch);
      }
      pch = strtok(NULL, ",");
      j++;
    }
    i++;
    current_data->next  = NULL;
        
  }
  current_data = root_data;
  model_at -> no_obs_sps = (j-1);/*minus 1 for time*/
  model_at -> no_obs = i;
   
  return(root_data);

}

