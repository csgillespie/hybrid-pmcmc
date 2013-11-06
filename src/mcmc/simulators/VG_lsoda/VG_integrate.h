extern void VG_integrate(double *y,double * t,double * tout,double * theta,int * odesize,
	void (* derivs) 
        (int * neq ,double *t, double *y, double *fout, double *theta),
	       double * reltol, double * abstol);


extern void lsoda_( void (* f) 
                      (int *,double *, double*,double*,double *),
		    const int *neq, double *y, double *t, double *tout, 
                    int *itol, double *rtol, double *atol, int *itask,
		    int *istate, int *iopt, double * rwork, int *lrw, 
		    int *iwork, int * liw,
		    void (*jac) 
                      (const int*, double*, double*,int*,int*,double*,int*), 
		    int *jt,double *params);

