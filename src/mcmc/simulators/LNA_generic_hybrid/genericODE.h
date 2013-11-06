// Functions for solving either LNA or 2MA equations.

void AllocReactionInfo(ReactionInfo *pI);
void DeAllocReactionInfo(ReactionInfo *pI);

void ODEInfoPrint(void *pOI);
void *pODEInfoAlloc(void *pRIv);
void ODEInfoDeAlloc(void *pOI);
void ODEInfoArrayPrint(int n, void **ppOI);
void **ppODEInfoArrayAlloc(int n_times, ReactionInfo *pRI);
void ODEInfoArrayDeAlloc(int n_times, void **ppOI);

void ODEInitialise(void *pOI,double t,double *pPostMean, double *pPostVar, 
		   double *peta);
void ODEIntegrateSimulateOne(double t, SDEInfo *pSI, void *pOI, double *pzero,
			     double *pxcurr, gsl_rng *r, double *pt_int, 
			     double *pt_no_int);
void ODEIntegrateOne(double t, SDEInfo *pSI, void *pOI);
void ODEExtract(void *pOI,double *pPriorMean, double *pPriorVar, double *peta,
		double *pt, double *ptmp); // ptmp is at least S*S

// Do not call this function yourself - it is called by the lsoda package
extern void lsoda_derivs(int *podedim, double *pt, double *py, double *pdydt, double *pparams);

