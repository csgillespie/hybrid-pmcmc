#include "../headers.h"
#include "mcmc_headers.h"

FILE* getFile(char filename[100]){

  FILE* out = fopen(filename, "w");
  if(NULL==out)
    {
      printf("Cannot open data file %s\n", filename);
      exit(1);
    }
  return(out);
}


void outputMCMCSettings(st_mcmc_settings *mcmc_settings)
{
  int thin, burn, parts;
  long iters;
    
  thin = mcmc_settings->thin;
  burn = mcmc_settings->burn;
  parts = mcmc_settings->parts;
  iters = mcmc_settings->iters;
    
  printf("#################################\n");
  printf("thin = %d\n", thin);
  printf("burn = %d\n", burn);
  printf("parts = %d\n", parts);
  printf("iters = %ld\n", iters);
  printf("#################################\n");
}

void outputPartSettings(st_part_at *accepted_part)
{
  printf("###########Part values############\n");
  printf("#params = %d\n", (int) accepted_part->params->size);
  printf("#sps = %d\n", (int) accepted_part->sps->size);
  printf("#res = %d\n", (int) accepted_part->res->size);
  printf("#################################\n");
}
    


void outputPartsSettings(st_parts_at *accepted_parts)
{
  printf("###########Parts values############\n");
  printf("#sps = %d, %d\n", (int) accepted_parts->sps->size1, (int) accepted_parts->sps->size2);
  printf("#res = %d, %d\n", (int) accepted_parts->res->size1, (int) accepted_parts->res->size2);
  printf("#################################\n");
}

void outputData(st_data *data_at, int no_obs_sps)
{
  int i;
  printf("#################################\n");
  printf("tstep = %f\n", data_at->tstep);
  for(i=0; i< no_obs_sps; i++)
    printf("%f ", data_at->obser[i]);
  printf("\n");
  printf("#################################\n");
}

void outputParameterVector(st_part_at *part)
{
  int  j;    
  for(j=0; j<part->params->size; j++)
    {
      printf("%3.3f ", VGET(part->params, j));
    }
  printf("%7.7f ", part->like);  
  printf("\n");
}

void outputParameters(int iteration, st_mcmc_settings *mcmc_settings, st_part_at *accepted_part)
{
  int  j;    
  FILE* out;
    
  out = mcmc_settings->output_file;
  /*    printf("%d %ld \n", accepted_parts->params->size2, mcmc_settings->iters);*/
  fprintf(out, "%d ", iteration);
  for(j=0; j< accepted_part->params->size; j++)
    {
      fprintf(out, "%7.10f ", VGET(accepted_part->params, j));
    }
  fprintf(out, "%7.7f ", accepted_part->like);  
  fprintf(out, "\n");    
  //fflush(out);
  /*      for(i=0; i<mcmc_settings->parts; i++){
          for(j=0; j< accepted_parts->sps->size2; j++)
          {
          fprintf(out, "%7.7f ", MGET(accepted_parts->sps, i, j));
          }*/
                
  /*}   */
}

void outputTimings(double t_diff, st_mcmc_settings *mcmc_settings)
{
  FILE* out;
  out = mcmc_settings->timing_file;
  fprintf(out, "%f\n", t_diff/1000.0);
  fflush(out);
}
