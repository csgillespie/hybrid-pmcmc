#include "mnormal.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include "model.h"
#ifdef __cplusplus
extern "C" {
#endif
  double logdens(double *initdata, double *data,int * nmols, double *tstart,
		 double *tend, double * thetas, int nthetas, int odesize, 
		 int  nmeans, double * relerr, double * abserr,
		 funderivs derivs,
		 spar ** mypar );
  double logdens2(double *initdata, double *data,int * nmols, double *tstart,
		  double *tend, double * thetas, int nthetas, int odesize, 
		  int nmeans,double * relerr, double * abserr,
		  funderivs derivs,
		 spar ** mypar );
  double logdens3(double *initdata, double *data,int * nmols, double *tstart,
		  double *tend, double * thetas, int nthetas, int odesize, 
		  int nmeans,double * relerr, double * abserr,
		  funderivs derivs,
		 spar ** mypar );
  double logdens4(spar ** mypar, double * thetas);

int readdata(const char * filename, spar ** mypar);
double llik(spar ** mypar, double *thetas);
gsl_vector * sample_points(double *initdata, double tstart,
	double tend, double * thetas, int nthetas, spar ** myparams, 
	 funderivs derivs
	);

SEXP calcdens(SEXP initdata, SEXP tstart,SEXP tend,
	      SEXP initmean, SEXP initvar,
	      SEXP thetas, SEXP s_odesize, 
	      SEXP relerr, SEXP abserr, 
	      SEXP ptrf);

#ifdef __cplusplus
}
#endif
