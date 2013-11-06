#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


// Scratch = workspace of S^2+S;
void SimulateTruncMVN(int S, double *pmean, double *pVar, 
		      double *pscratch, double *pZ, gsl_rng *r) {
  int i,j;
  double *pVarroot=pscratch;
  double *piid=pVarroot+S*S;

  matrix_square_root(S, pVar, pVarroot);

  for (i=0;i<S;i++) {
    piid[i]=gsl_ran_gaussian(r,1.0);
    pZ[i]=pmean[i];
  }

  for (i=0;i<S;i++) {
    for (j=0;j<S;j++) {
      pZ[i] += pVarroot[i*S+j]*piid[j];
      if (pZ[i]<0) {
	pZ[i]=0;
      }
    }
  }
}
