#include <stdio.h>
#include "lik.h"
#include "integrate.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
//#define GSLINTEGRATE
//#define DEBUG

void dointegration(double *y, double *t, double * t1, double * thetas,
	int * odesize, void*  derivs, double *rerr, double *aerr) {
#ifdef GSLINTEGRATE
    /* start the GSL numerical integrator */
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;

    //const gsl_odeiv_step_type * T = gsl_odeiv_step_bsimp;
    const gsl_odeiv_step_type * T1 = gsl_odeiv_step_rkf45;
    //const gsl_odeiv_step_type * T = gsl_odeiv_step_rkck;
    //const gsl_odeiv_step_type * T = gsl_odeiv_step_gear2;
    //const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
    
    gsl_odeiv_step * s;
    if(myparams->jac==NULL)
      {
	s = gsl_odeiv_step_alloc (T1, odesize);
      }
    else
      {
	s = gsl_odeiv_step_alloc (T, odesize);
      }
    gsl_odeiv_control * c = gsl_odeiv_control_y_new ( aerr, rerr);
    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (odesize);
    gsl_odeiv_system sys = {derivs,myparams->jac , odesize, thetas};


    while (t < t1)
    {
	int status = gsl_odeiv_evolve_apply (e, c, s,
		&sys,
		&t, t1,
		 rerr, y);

	if (status != GSL_SUCCESS) break;
#ifdef DEBUG
	printf("t=%f ",t);
	for(i=0;i<odesize;i++) printf("%4.9f ",y[i]);
 	printf("\n");
#endif
    }
    /*
     * Initialize vectors needed for the multivariate normal log-density
     */
    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free(s);

#else
    integrate(y, t, t1,thetas, odesize,derivs, rerr, aerr);
#endif
}



double llik(spar ** mypar, double *thetas) 
{
    spar * myparams=*mypar;
    int minrows = myparams->nrows -1;
    int mincols = myparams->ncols -1;
    register int i;
    double sum=0;
    switch(myparams->method){
    case 1:  /* Concetrations */
      for(i=0;i<minrows;i++)
	{
	  sum=logdens( (myparams->dataset) +1+ (i* (myparams->ncols)), /* initdata */
		       (myparams->dataset) +1+ ( (i+1) * (myparams->ncols)), /* data */
		       & myparams->nmols,
		       (myparams->dataset) + (i* (myparams->ncols)),
		       (myparams->dataset) + ((i+1)* (myparams->ncols)),thetas,
		       myparams->sysdim,myparams->odesize,mincols,
		       &myparams->relerr, &myparams->abserr, myparams->f,mypar)+sum;
	};
      break;
    /* Molecules */
    case 0:
      for(i=0;i<minrows;i++)
	{
	  
	  sum=logdens2( (myparams->dataset) +1+ (i* (myparams->ncols)), /* initdata */
			(myparams->dataset) +1+ ( (i+1) * (myparams->ncols)), /* data */
			& myparams->nmols,
			(myparams->dataset) + (i* (myparams->ncols)),
			(myparams->dataset) + ((i+1)* (myparams->ncols)),thetas,
			myparams->sysdim,myparams->odesize,mincols,
			&myparams->relerr, &myparams->abserr, myparams->f,mypar)+sum;
	};
      break;
    /* Concentrations and non zero M_t  */
    case 3: 
      for(i=0;i<minrows;i++)
	{
	  
	  sum=logdens3( (myparams->dataset) +1+ (i* (myparams->ncols)),       /* Initial Observations */
			(myparams->dataset) +1+ ( (i+1) * (myparams->ncols)), /* Subsequent Observations */
			& myparams->nmols,                                     /* System Size */
			(myparams->dataset) + (i* (myparams->ncols)),         /* Starting Time */
			(myparams->dataset) + ((i+1)* (myparams->ncols)),     /* Ending Time */
			thetas,            /* Thetas array */
			myparams->sysdim,  /* Number of thetas */
			myparams->odesize, /* Size of Odes*/
			mincols,           /* Dimension of the ODEs of Means (Y_t) */
			&myparams->relerr, /* Integrator's relative error */
			&myparams->abserr, /* Integrator's absolute error */
			myparams->f,       /* Pointer to the derivatives function */
			mypar              /* Pointer to the model Structure */
			)+sum;
	};
      break;
    /* Non restarting method */
    case 4:
      sum=logdens4(mypar,thetas);
      break;
    }
    if(isnan(sum)) {
	fprintf(stderr,"Log-likelihood evaluation error - try a different initial value\n");
	printf("Log-likelihood evaluation error - try a different initial value\n");
	printf("Thetas:\n");
	for(i=0;i<myparams->x->size;i++) printf("  %f\n",thetas[i]);
	exit(10);
    }
    return sum;
}

/************************************************ 
   transition log-density based on Concentrations 
**************************************************/
double logdens(double *initdata, double *data,int * nmols, double *tstart,
	       double *tend, double * thetas, int nthetas, int odesize, 
	       int nmeans, double  * relerr, double  * abserr, 
               funderivs derivs,
	       spar ** mypar
	       )
{
    int i,j,sysind;
    double y[odesize],rerr,aerr,ntot;
    ntot = *nmols;
    for(i=0; i< nmeans; i++) { y[i]=initdata[i]/ ntot ;}
    for(i=nmeans;i<odesize;i++) { y[i]=0; }
    double t = *tstart , t1 = *tend;
    rerr=*relerr;aerr=*abserr;
    dointegration(y,&t,&t1,thetas,&odesize,(void *) derivs,&rerr,&aerr);

    /* turn y to mean vector - var-cov matrix */
    
    sysind=0;
    gsl_vector * meanv  = gsl_vector_alloc (nmeans);
    gsl_vector * datav  = gsl_vector_alloc (nmeans);
    for (i=0;i<nmeans;i++)
    {
	gsl_vector_set(meanv,i,y[i]* (*nmols) );
	gsl_vector_set(datav,i,data[i]);
	sysind++;
    }


    gsl_matrix * covm = gsl_matrix_alloc (nmeans, nmeans);

    for (i=0;i<nmeans;i++)
    {
	for (j=i;j<nmeans;j++)
	{
	    gsl_matrix_set (covm, i, j, y[sysind]* (*nmols) );
	    gsl_matrix_set (covm, j, i, y[sysind]* (*nmols) );
	    sysind++;
	}
    }

    double dens = ldmvnorm( nmeans, datav, meanv, covm);
#ifdef DEBUG
    for (i=0;i<nmeans;i++)
    {
	printf("%f ",gsl_vector_get(meanv,i));
    }
    printf("\n");
    for (i=0;i<nmeans;i++)
    {
	printf("%f ",gsl_vector_get(datav,i));
    }
    printf("\n");
    for (i=0;i<nmeans;i++)
    {
	for (j=0;j<nmeans;j++)
	{
	    printf("%f ",gsl_matrix_get(covm,i,j));
	}
	printf("\n");
    }
#endif

    gsl_vector_free(meanv);
    gsl_vector_free(datav);
    gsl_matrix_free (covm);
    return dens;
}

/***********************************************
   transition log-density based on Concentrations
   Assuming Non-Zero M-elements
**************************************************/
double logdens3(double *initdata, double *data,int * nmols, double *tstart,
	       double *tend, double * thetas, int nthetas, int odesize, 
	       int nmeans, double  * relerr, double  * abserr, 
               funderivs derivs,
	       spar ** mypar
	       )
{
    spar * myparams =  *mypar;
    int i,j,sysind;
    double y[odesize],rerr,aerr;

    double t , t1;
    double ntot= *nmols;
    double sntot= sqrt(ntot);

    t=myparams->dataset[0]; t1=*tstart;

    rerr=*relerr;
    aerr=*abserr;
    for(i=0; i< nmeans; i++) { y[i]=myparams->dataset[i+1]/ ntot ;}
    for(i=nmeans;i<odesize;i++) { y[i]=0; }

    if(t1-t>1e-5)
    {
	for(i=0; i< nmeans; i++) { y[i]=myparams->dataset[i+1]/ ntot ;}
	for(i=nmeans;i<odesize;i++) { y[i]=0; }
	dointegration(y,&t,&t1,thetas,&odesize,(void *) derivs,&rerr,&aerr);

	/*Reset Integration parameters */
	rerr=*relerr;
	aerr=*abserr;
	for(i=nmeans;i<odesize;i++) { y[i]=0; }
	for(i=0;i<nmeans;i++) 
	{ 
	    y[i+odesize-nmeans]=initdata[i]/sntot - y[i]*sntot; 
	}
    }

    t = *tstart; t1 = *tend;

    dointegration(y,&t,&t1,thetas,&odesize,(void *) derivs,&rerr,&aerr);

    /* turn y to mean vector - var-cov matrix */
    
    sysind=0;
    gsl_vector * meanv  = gsl_vector_alloc (nmeans);
    gsl_vector * datav  = gsl_vector_alloc (nmeans);
    for (i=0;i<nmeans;i++)
    {
        gsl_vector_set(meanv,i,y[i] * ntot +y[i+14] * sntot );
	//gsl_vector_set(meanv,i,y[i] * (*nmols) );
	gsl_vector_set(datav,i,data[i]);
	sysind++;
    }


    gsl_matrix * covm = gsl_matrix_alloc (nmeans, nmeans);

    for (i=0;i<nmeans;i++)
    {
	for (j=i;j<nmeans;j++)
	{
	    gsl_matrix_set(covm, i, j, y[sysind]* ntot );
	    gsl_matrix_set(covm, j, i, y[sysind]* ntot );
	    sysind++;
	}
    }

#ifdef DEBUG
    for (i=0;i<nmeans;i++)
    {
	printf("%f ",gsl_vector_get(meanv,i));
    }
    printf("\n");
    for (i=0;i<nmeans;i++)
    {
	printf("%f ",gsl_vector_get(datav,i));
    }
    printf("\n");
    for (i=0;i<nmeans;i++)
    {
	for (j=0;j<nmeans;j++)
	{
	    printf("%f ",gsl_matrix_get(covm,i,j));
	}
	printf("\n");
    }
#endif

    double dens = ldmvnorm( nmeans, datav, meanv, covm);
    gsl_vector_free(meanv);
    gsl_vector_free(datav);
    gsl_matrix_free (covm);
    return dens;
}


/*
 * transition log-density using number of molecules
 * */
double logdens2(double *initdata, double *data,int * nmols, double *tstart,
		double *tend, double * thetas, int nthetas, int odesize, 
		int nmeans, double * relerr, double * abserr, 
                funderivs derivs,
		spar ** mypar
		)
{
    int i,j,sysind;
    double y[odesize],rerr,aerr;
    double t = *tstart , t1 = *tend;

    for(i=0;i<nmeans;i++) { y[i]=initdata[i];}
    for(i=nmeans;i<odesize;i++) { y[i]=0; }

    rerr=*relerr;aerr=*abserr;
    dointegration(y,&t,&t1,thetas,&odesize,(void *) derivs, &rerr, &aerr);
    
    sysind=0;
    gsl_vector * meanv  = gsl_vector_alloc (nmeans);
    gsl_vector * datav  = gsl_vector_alloc (nmeans);
    for (i=0;i<nmeans;i++)
    {
	gsl_vector_set(meanv,i,y[i]);
	gsl_vector_set(datav,i,data[i]);
	sysind++;
    }

    gsl_matrix * covm = gsl_matrix_alloc(nmeans, nmeans);

    for (i=0;i<nmeans;i++)
    {
	for (j=i;j<nmeans;j++)
	{
	    gsl_matrix_set (covm, i, j, y[sysind]);
	    gsl_matrix_set (covm, j, i, y[sysind]);
	    sysind++;
	}
    }

    #ifdef DEBUG
    for (i=0;i<nmeans;i++)
    {
	printf("%f ",gsl_vector_get(meanv,i));
    }
    printf("\n");
    for (i=0;i<nmeans;i++)
    {
	printf("%f ",gsl_vector_get(datav,i));
    }
    printf("\n");

    for (i=0;i<nmeans;i++)
    {
	for (j=0;j<nmeans;j++)
	{
	    printf("%f ",gsl_matrix_get(covm,i,j));
	}
	printf("\n");
    }
    #endif
    double dens = ldmvnorm( nmeans, datav, meanv, covm);

    gsl_vector_free(meanv);
    gsl_vector_free(datav);
    gsl_matrix_free (covm);
    return dens;
}
gsl_vector * sample_points(double *initdata, double tstart,
	double tend, double * thetas, int nthetas, spar ** mypar, 
	 funderivs derivs
	)
{
    spar * myparams = (spar *) *mypar;
    int nmols=myparams->nmols;
    int odesize=myparams->odesize;
    int nmeans=myparams->odesize-myparams->sysdim;
    double relerr=myparams->relerr;
    double abserr=myparams->abserr;
    int i,j,sysind;
    double y[odesize];
    for(i=0;i<nmeans;i++) { y[i]=initdata[i]/nmols;}
    for(i=nmeans;i<odesize;i++) { y[i]=0; }

    /* start the numerical integrator */
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
    //const gsl_odeiv_step_type * T = gsl_odeiv_step_rkck;
    //const gsl_odeiv_step_type * T = gsl_odeiv_step_gear2;
    //const gsl_odeiv_step_type * T = gsl_odeiv_step_rk8pd;
    gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, odesize);
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (relerr, abserr);
    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (odesize);


    gsl_odeiv_system sys = {derivs, 0, odesize, thetas};
 
    double t = tstart , t1 = tend;

    while (t < t1)
    {
	int status = gsl_odeiv_evolve_apply (e, c, s,
		&sys, 
		&t, t1,
		&relerr, y);

	if (status != GSL_SUCCESS)
	    break;

    }

    /*
     * Initialize vectors needed for the multivariate normal log-density
     */
    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);
    /* turn y to mean vector - var-cov matrix */
    
    sysind=0;
    gsl_vector * meanv  = gsl_vector_alloc (nmeans);
    for (i=0;i<nmeans;i++)
    {
	gsl_vector_set(meanv,i,y[i]*nmols);
	sysind++;
    }


    gsl_matrix * covm = gsl_matrix_alloc (nmeans, nmeans);

    for (i=0;i<nmeans;i++)
    {
	for (j=i;j<nmeans;j++)
	{
	    gsl_matrix_set (covm, i, j, y[sysind]*nmols);
	    gsl_matrix_set (covm, j, i, y[sysind]*nmols);
	    sysind++;
	}
    }

    gsl_vector * gpoints = gsl_vector_alloc(nmeans);
    rmvnorm( r, nmeans, meanv, covm,gpoints);
    gsl_vector_free(meanv);
    gsl_matrix_free (covm);
    return gpoints;
}

/*********************************************** 
   Non-restarting method 
************************************************/
double logdens4(spar ** mypar, double * thetas)
{

    spar * myparams =  *mypar;
    int i,j,k,sysind,mindex;
    int nmeans=  myparams->ncols-1;
    double y[myparams->odesize],rerr,aerr;
    double t , t1;
    double ntot= myparams->nmols;
    double sntot= sqrt(ntot);
    double sum=0;

    /*
     *y vector has the following composition:
      y = Y_0,...,Y_nmeans,S_[1 1],...,S_[nmeans nmeans],
          M_0,...,M_nmeans

      Y: ODE corresponding to the macroscopic approximation
      S: upper diagonal of covariance matrix
      M: mean of linear SDE
     */

    /* Points to the position of the M vector*/
    mindex=myparams->odesize-nmeans;


    /* GSL variables Declarations */
    gsl_vector * meanv  = gsl_vector_alloc (nmeans);
    gsl_vector * datav  = gsl_vector_alloc (nmeans);
    gsl_matrix * covm = gsl_matrix_alloc (nmeans, nmeans);
    /* End of GSL Declarations */

    t=myparams->dataset[0]; t1=myparams->dataset[ myparams->ncols];

    rerr= myparams->relerr;
    aerr= myparams->abserr;

    /* Initialize considering M_t0 = (0,0, ... ,0)*/
    for(i=0; i< nmeans; i++) 
    { 
	y[i]=myparams->dataset[i+1]/ ntot;
    }
    /* Initialize M_t and S_t */
    for(i=nmeans;i<myparams->odesize;i++) { y[i]=0; }

    /*Solve the ODEs*/
    dointegration(y,&t,&t1,thetas, &myparams->odesize,myparams->f,&rerr,&aerr);
    

    sysind=0;
    for (i=0;i<nmeans;i++)
    {
	gsl_vector_set(meanv,i,y[i] * ntot +y[i+mindex] * sntot );
	gsl_vector_set(datav,i,myparams->dataset[ myparams->ncols + i+1]);
	sysind++;
    }
    
    for (i=0;i<nmeans;i++)
    {
	for (j=i;j<nmeans;j++)
	{
	    gsl_matrix_set(covm, i, j, y[sysind]* ntot );
	    gsl_matrix_set(covm, j, i, y[sysind]* ntot );
	    sysind++;
	}
    }
    sum = ldmvnorm( nmeans, datav, meanv, covm);

    for(k=1;k< (myparams->nrows - 1);k++)
    {
	/*Initialize S_t */
	for(i=nmeans;i< mindex;i++) { y[i]=0; }

	/*Initialize M_t*/
	for(i=0;i<nmeans;i++) 
	{ 
	    y[i+mindex]=
		myparams->dataset[k*(myparams->ncols)+i+1]/sntot - y[i]*sntot; 
	    
	}
	t=myparams->dataset[k* (myparams->ncols)];
	t1=myparams->dataset[(k+1)* (myparams->ncols)];


	/*Reset Integration parameters */
	rerr= myparams->relerr;
	aerr= myparams->abserr;

	/*Solve forward the ODEs*/
	dointegration(y,&t,&t1,thetas,& myparams->odesize,
		(void *) myparams->f,&rerr,&aerr);

	/* Calculate the transition probability*/
	sysind=0;
	for (i=0;i<nmeans;i++)
	{
	    gsl_vector_set(meanv,i,y[i] * ntot +y[i+mindex] * sntot );
	    gsl_vector_set(datav,i,myparams->dataset[(k+1)*(myparams->ncols)+i+1]);
	    sysind++;
	}

	for (i=0;i<nmeans;i++)
	{
	    for (j=i;j<nmeans;j++)
	    {
		gsl_matrix_set(covm, i, j, y[sysind]* ntot );
		gsl_matrix_set(covm, j, i, y[sysind]* ntot );
		sysind++;
	    }
	}

#ifdef DEBUG
    printf("Means ");
    for (i=0;i<nmeans;i++)
    {
	printf("%f ",gsl_vector_get(meanv,i));
    }
    printf("\n--\n");
    printf("Data ");
    for (i=0;i<nmeans;i++)
    {
	printf("%f ",gsl_vector_get(datav,i));
    }
    printf("\n--\n");
    printf("Cov\n ");
    for (i=0;i<nmeans;i++)
    {
	for (j=0;j<nmeans;j++)
	{
	    printf("%f ",gsl_matrix_get(covm,i,j));
	}
	printf("\n");
    }
	printf("--\n");
#endif
	sum += ldmvnorm( nmeans, datav, meanv, covm);
	
    }

#ifdef DEBUG
    printf("Means ");
    for (i=0;i<nmeans;i++)
    {
	printf("%f ",gsl_vector_get(meanv,i));
    }
    printf("\n");
    printf("Data ");
    for (i=0;i<nmeans;i++)
    {
	printf("%f ",gsl_vector_get(datav,i));
    }
    printf("\n");
    printf("Cov\n");
    for (i=0;i<nmeans;i++)
    {
	for (j=0;j<nmeans;j++)
	{
	    printf("%f ",gsl_matrix_get(covm,i,j));
	}
	printf("\n");
    }
#endif

    gsl_vector_free(meanv);
    gsl_vector_free(datav);
    gsl_matrix_free (covm);
    return sum;
}


/*************************************************************
   Calculates the transition density and its parameters 
   based on Concentrations 
**************************************************************/
SEXP calcdens(SEXP initdata,
	      SEXP tstart,SEXP tend,
	      SEXP initmean, SEXP initvar,
	      SEXP thetas, SEXP sodesize, 
	      SEXP  relerr, SEXP abserr, 
	      SEXP ptrf)
{
    int i,j,sysind,k,z;
    int odesize = INTEGER(sodesize)[0];
    int nmeans = length(initdata);
    int lentimes = length(tend);
    double y[odesize],rerr,aerr;
    funderivs derivs = (funderivs ) R_ExternalPtrAddr(ptrf);
    double *p_thetas = NUMERIC_POINTER(thetas);
    /* Initialize the vector to be integrated */
    for(i=0; i< nmeans; i++) 
    { 
      y[i]=REAL(initdata)[i]; 
    //  Rprintf("%f\n",y[i]);
    }
    j=odesize-nmeans;
    z=0;
    for(i=nmeans;i<j;i++) 
    { 
      y[i]=REAL(initvar)[z]; 
      //Rprintf("%f\n",y[i]);
      z++;
    }
    z=0;
    for(i=j;i<odesize;i++)
    { 
      y[i]=REAL(initmean)[z]; 
      //Rprintf("%f\n",y[i]);
      z++;
    }
    
    double *p_tend=NUMERIC_POINTER(tend);
    double t= REAL(tstart)[0],t1;
    rerr= REAL(relerr)[0];
    aerr=REAL(abserr)[0];

    /* Used as the return elements of the list */
    SEXP meanv,odev,covm,time;

    //Point to the named elements of the timepoint list
    double *p_time, *p_meanv,*p_odev, *p_covm; 

    SEXP list_names,list,ans;
    char *names[4]={"Time","ODE","MEAN","VAR"};

    /*Pass the names in a char vector */
    PROTECT(list_names = allocVector(STRSXP, 4));
    for(i = 0; i < 4; i++)
      SET_STRING_ELT(list_names, i,  mkChar(names[i]));
    
    /*Create an empty list for each timepoint*/
    PROTECT(ans = allocVector(VECSXP, lentimes));


    /*-----------------------------------------------------*/
    /* The loop starts */
    for(k=0;k<lentimes;k++){

      PROTECT(time = NEW_NUMERIC(1));
      p_time = NUMERIC_POINTER(time);
      PROTECT(meanv = NEW_NUMERIC(nmeans));
      p_meanv = NUMERIC_POINTER(meanv);      
      PROTECT(odev = NEW_NUMERIC(nmeans));
      p_odev = NUMERIC_POINTER(odev);
      PROTECT(covm = allocMatrix(REALSXP,nmeans,nmeans));
      p_covm = NUMERIC_POINTER(covm);


      double t1 = p_tend[k];
      *p_time = t1;
      if(k==0)
	integrate(y,&t,&t1,p_thetas,&odesize, derivs,&rerr,&aerr);
      else //Resumes the integration
	{
	  t=REAL(tend)[k-1];
	  rerr= REAL(relerr)[0];
	  aerr=REAL(abserr)[0];
	  integrate(y,&t,&t1,p_thetas,&odesize, derivs,&rerr,&aerr);
	  //cont_integr(y,&t,&t1,p_thetas,&odesize,(void *) derivs,&rerr,&aerr);
	}
      
      for (i=0;i<nmeans;i++)
	{
	  p_odev[i] = y[i];
	  p_meanv[i] = y[i+odesize-nmeans];
	}
      sysind=nmeans;
      
      for (i=0;i<nmeans;i++)
	{
	  for (j=i;j<nmeans;j++)
	    {
	      p_covm[i + nmeans * j] = y[sysind];
	      p_covm[i * nmeans + j] = y[sysind];
	      sysind++;
	    }
	}
                     
      PROTECT(list = allocVector(VECSXP, 4)); //Alocate list's elements
      SET_VECTOR_ELT(list, 0, time);
      SET_VECTOR_ELT(list, 1, odev);
      SET_VECTOR_ELT(list, 2, meanv);
      SET_VECTOR_ELT(list, 3, covm);
      setAttrib(list, R_NamesSymbol, list_names); //Set list's names
      SET_VECTOR_ELT(ans,k,list);
      UNPROTECT(5);
    }
    UNPROTECT(2);
    return(ans);
}
