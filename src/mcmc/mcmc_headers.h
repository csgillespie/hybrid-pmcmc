

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>



#define VSET gsl_vector_set
#define VGET gsl_vector_get
#define MSET gsl_matrix_set
#define MGET gsl_matrix_get



typedef struct
{
  int thin;
  int burn;
  long iters;
  int parts;
  FILE *output_file;
  FILE *timing_file;/*In millisecs*/
} st_mcmc_settings;



/* General look up table about the model */
typedef struct
{
  int no_obs; //# of observations
  int no_sps; //# of species
  int no_obs_sps; //# of observed species
  int no_params; //# of parameters
} st_model_at;


/* The data structure is a linked list. In theory this 
   makes it easier to have different inter-observation times */
typedef struct dataStruct
{
  double *obser;
  double tstep;
  struct dataStruct *next;
} st_data;


/* 
 * A particular particle
 * Under PMCMC params is the same for all particles
 */
typedef struct
{
  gsl_vector *params;
  gsl_vector *sps;
  gsl_vector *res;/*used only in the hybrid simulator*/
  double like;
} st_part_at;


/*Particle matrix */
typedef struct
{
  gsl_matrix *sps;
  gsl_matrix *res;/*used only in the hybrid simulator*/
} st_parts_at;



/* always_accept means that every proposal is accepted.
 * Always accept occurs when all parameters are fixed.
 * Used when estimating the number of particles when 
 * all parameters are fixed 
 */
typedef struct
{
  gsl_matrix *tuning; /*Cholsky decomp */
  gsl_vector *fixed;/* Which parameters vary */
  gsl_vector *z; /*used to store N(0, 1)*/
  int always_accept;
} st_mcmc_update;

