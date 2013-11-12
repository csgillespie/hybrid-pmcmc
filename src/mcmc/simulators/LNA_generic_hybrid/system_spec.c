#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#include "struct_generic.h"
#include "../Generic_c/LinearAlgebra/cslinalg.h"


// User changeable variables:
// csns:      number of (free - unconstrained) species
// csnr:      number of reactions
// nits:      # iterations of RWM
// seed:      random seed
// psystem:   descriptive string about the system
// pobsfname: file name for obs
// poutfname: file name for mcmc output

// User changeable functions:
// get_A:    net effects matrix - one row for each reaction
// get_h:    reaction rates vector from rate constants and species amounts
// get_dhdx: matrix of derivative of reaction rates vector wrt each species
// get_d2hdxmdxn: hessian of the hcpt component of h - 2MA only
// get_fixed_parms: any fixed parameters (i.e. not rate constants or err vars)
//                  or returns zero if there are none
// constrain: correction if ODE or Kalman takes species outside allowed range
//            is called by the main code before get_h etc.
// get_prior_for_X0: apriori species at t=0 are MVN; supply the mean and var.
// initialise_rates: starting values for vector (returns log(c)) 
// initialise_obsvar: starting (or fixed) values for the matrix
// allowable_rates: are all rate values in permitted range?  0/1 - RWM
// allowable_obsvar: are all obs variance values in permitted range?  0/1 RWM
// log_prior_rates: given the rate vec return the log prior for LOG of the 
//                  parameter. RWM.
// log_prior_obsvar: given the obs var mx return the log prior RWM
// get_jumpvar: variance mx for log of the parameter RWM


// Other functions - don't change these!
// free_fixed_params;
// get_SR;
// get_nits_and_seed
// pget_system
// pget_obsfname
// pget_outfname

// INSTRUCTION: choose the system by commenting out the unnecessary #defines 

//#define GGA
#define GGB

static void ForceNonNeg(double *pv, int n) {
  int i;
  for (i=0;i<n;i++) {
    pv[i]*=(pv[i]>0);
  }
}

static int hybrid_use_sk=0; // 0= use new decision method, 1=S&K
static double hybrid_hstar=4; // S&K

static double hybrid_N=60; // CS and S&K (this is N/epsilon in the paper)

// change this to alter the frequency of checks on which reactions are 
// fast and which are slow
static double hybrid_deltat=0.1; 
static double hybrid_epsilon=0.25; 



#ifdef GGA
// 2 species and 5 reactions
static int csns=2;
static int csnr=5;
static int nits=1000;
static int seed=1232123; 
static char psystem[]="GG 2 species\n";
//static char pobsfname[]="../../Obs/GG100.obs";
//static poutfname[]="GG100x100000.mcmc";
static char pobsfname[]="../../Obs/SIR.obs";
static char poutfname[]="SIRx1000.mcmc";

//char poutfname[]="blob.mcmc";

//0  0->X1
//1  0->X2
//2  X1->0
//3  X2->0
//4  X1+X2->2X2

// SIR is special case with X1=S, X2=I and just reactions 3 and 4.

void get_A(int *pA, int S, int R) {
  int pA0[]={1,0};
  int pA1[]={0,1};
  int pA2[]={-1,0};
  int pA3[]={0,-1};
  int pA4[]={-1,1};

  int i;

  assert(S==csns);
  assert(R==csnr);

  for (i=0;i<csns;i++) {
    pA[0*csns+i]=pA0[i];
    pA[1*csns+i]=pA1[i];
    pA[2*csns+i]=pA2[i];
    pA[3*csns+i]=pA3[i];
    pA[4*csns+i]=pA4[i];
  }
}

void constrain(double *px) {

  if (px[0]<0) {
    if (px[0]<-0.001) {
      printf("x0 was %f\n",*px);
    }
    px[0]=0;
  }
  if (px[1]<0) {
    if (px[1]<-0.001) {
      printf("x1 was %f\n",px[1]);
    }
    px[1]=0;
  }
}


void get_h(double *ph, double *px, double *pc, void *pparms) {
  ph[0]=pc[0];
  ph[1]=pc[1];
  ph[2]=pc[2]*px[0];
  ph[3]=pc[3]*px[1];
  ph[4]=pc[4]*px[0]*px[1];
  ForceNonNeg(ph,csnr);
}

void get_dhdx(double *pdhdx, double *px, double *pc, void *pparms) {
  pdhdx[0*csns+0]=0; pdhdx[0*csns+1]=0; 
  pdhdx[1*csns+0]=0; pdhdx[1*csns+1]=0; 
  pdhdx[2*csns+0]=pc[2]; pdhdx[2*csns+1]=0; 
  pdhdx[3*csns+0]=0; pdhdx[3*csns+1]=pc[3]; 
  pdhdx[4*csns+0]=pc[4]*px[1]; 
  pdhdx[4*csns+1]=pc[4]*px[0];
}

// Linear term in Delta x in Taylor expansion of h(x+ \Delta x)
void get_b(double *pb, double *px, double *pc) {
  pb[0]=pc[2]+pc[4]*px[1];
  pb[1]=pc[3]+pc[4]*px[0];
}

void get_d2hdxmdxn(double *pd2hdxmdxn, double *px, double *pc, void *pparms,int hcpt) {
  vec_set(pd2hdxmdxn,4,0.0);

  switch(hcpt) {
  case 4: {
    pd2hdxmdxn[0*2+1]=pc[4];
    pd2hdxmdxn[1*2+0]=pc[4];
    break;
  }
  }
}

void* get_fixed_parms() {
  return 0;
}

void get_prior_for_X0(double * pmean, double *pvar) {
  int i,j;

  for (i=0;i<csns;i++) {
    pmean[i]=100;
    for (j=0;j<csns;j++) {
      if (i==j) {
	pvar[i*csns+j]=2500;
      }
      else {
	pvar[i*csns+j]=0.0;
      }
    }
  }
}

void initialise_rates(double *pc) {
  // GGS 2011 paper
  pc[0]=10; pc[1]=0.1; pc[2]=0.1; pc[3]=0.7; pc[4]=0.008;
  //  pc[0]=100; pc[1]=1.0;
}

void initialise_obsvar(int Sobs, double *pobsvar) {
  int i,j;
  for (i=0;i<Sobs;i++) {
     for (j=0;j<Sobs;j++) {
      if (i==j) {
	pobsvar[i*Sobs+j]=0.000001;
      }
      else {
	pobsvar[i*Sobs+j]=0.0;
      }
    }
  }
}

// Uniform on the log scale between -8 and 8
int allowable_rates(double *pc){
  int i, allowable=1;
  for (i=0;i<csnr;i++) {
    if ((pc[i]<exp(-8)) || (pc[i]>exp(8))) {
      allowable=0;
    }
  } 

  return(allowable);
}
int allowable_obsvar(int Sobs, double *pobsvar) {
  int i, allowable=1;
  for (i=0;i<Sobs;i++) {
    if ((pobsvar[i*Sobs+i]<-exp(20)) || (pobsvar[i*Sobs]>exp(3))) {
      allowable=0;
    }
  }
  return(allowable);
}

double log_prior_rates(double *pc) {
  return 0;
}

double log_prior_obsvar(int Sobs, double *pobsvar) {
  double lp=0;
  int i;

  for (i=0;i<Sobs;i++) {
    lp -= log(1+pobsvar[i*Sobs+i]*pobsvar[i*Sobs+i]); // Cauchy
  }

  return(lp); 
}

void get_jumpvar(int Sobs, double *pjv) {
  int i,j;

  // Zero first
  for (i=0;i<csnr+Sobs;i++) {
    for (j=0;j<csnr+Sobs;j++) {
	pjv[i*(Sobs+csnr)+j]=0.0;
    }
  }
  printf("%d %d\n",csnr,Sobs);
  for (i=0;i<csnr;i++) {
    pjv[i*(csnr+Sobs)+i]=0.01;
  }
}



#endif


#ifdef GGB // only difference is that final reaction produces 21 X2
// 2 species and 5 reactions
static int csns=2;
static int csnr=5;
static int nits=1000; // irrelevant as no mcmc here
static int seed=1232123; 
static char psystem[]="GGfriendly 2 species\n";
//static char pobsfname[]="../../Obs/GG100.obs";
//static poutfname[]="GG100x100000.mcmc";
static char pobsfname[]="../../Obs/SIR.obs"; // irrelevant
static char poutfname[]="SIRx1000.mcmc"; // ditto

//char poutfname[]="blob.mcmc";

//0  0->X1
//1  0->X2
//2  X1->0
//3  X2->0
//4  X1+X2->21X2

// SIR is special case with X1=S, X2=I and just reactions 3 and 4.

void get_A(int *pA, int S, int R) {
  int pA0[]={1,0};
  int pA1[]={0,1};
  int pA2[]={-1,0};
  int pA3[]={0,-1};
  int pA4[]={-1,20};

  int i;

  assert(S==csns);
  assert(R==csnr);

  for (i=0;i<csns;i++) {
    pA[0*csns+i]=pA0[i];
    pA[1*csns+i]=pA1[i];
    pA[2*csns+i]=pA2[i];
    pA[3*csns+i]=pA3[i];
    pA[4*csns+i]=pA4[i];
  }
}

void constrain(double *px) {

  if (px[0]<0) {
    if (px[0]<-0.001) {
      printf("x0 was %f\n",*px);
    }
    px[0]=0;
  }
  if (px[1]<0) {
    if (px[1]<-0.001) {
      printf("x1 was %f\n",px[1]);
    }
    px[1]=0;
  }
}


void get_h(double *ph, double *px, double *pc, void *pparms) {
  ph[0]=pc[0];
  ph[1]=pc[1];
  ph[2]=pc[2]*px[0];
  ph[3]=pc[3]*px[1];
  ph[4]=pc[4]*px[0]*px[1];
  ForceNonNeg(ph,csnr);
}

void get_dhdx(double *pdhdx, double *px, double *pc, void *pparms) {
  pdhdx[0*csns+0]=0; pdhdx[0*csns+1]=0; 
  pdhdx[1*csns+0]=0; pdhdx[1*csns+1]=0; 
  pdhdx[2*csns+0]=pc[2]; pdhdx[2*csns+1]=0; 
  pdhdx[3*csns+0]=0; pdhdx[3*csns+1]=pc[3]; 
  pdhdx[4*csns+0]=pc[4]*px[1]; 
  pdhdx[4*csns+1]=pc[4]*px[0];
}

// Linear term in Delta x in Taylor expansion of h(x+ \Delta x)
void get_b(double *pb, double *px, double *pc) {
  pb[0]=pc[2]+pc[4]*px[1];
  pb[1]=pc[3]+pc[4]*px[0];
}

void get_d2hdxmdxn(double *pd2hdxmdxn, double *px, double *pc, void *pparms,int hcpt) {
  vec_set(pd2hdxmdxn,4,0.0);

  switch(hcpt) {
  case 4: {
    pd2hdxmdxn[0*2+1]=pc[4];
    pd2hdxmdxn[1*2+0]=pc[4];
    break;
  }
  }
}

void* get_fixed_parms() {
  return 0;
}

void get_prior_for_X0(double * pmean, double *pvar) {
  int i,j;

  for (i=0;i<csns;i++) {
    pmean[i]=100;
    for (j=0;j<csns;j++) {
      if (i==j) {
	pvar[i*csns+j]=2500;
      }
      else {
	pvar[i*csns+j]=0.0;
      }
    }
  }
}

void initialise_rates(double *pc) {
  // GGS 2011 paper
  pc[0]=10; pc[1]=0.1; pc[2]=0.1; pc[3]=0.7; pc[4]=0.008;
  //  pc[0]=100; pc[1]=1.0;
}

void initialise_obsvar(int Sobs, double *pobsvar) {
  int i,j;
  for (i=0;i<Sobs;i++) {
     for (j=0;j<Sobs;j++) {
      if (i==j) {
	pobsvar[i*Sobs+j]=0.000001;
      }
      else {
	pobsvar[i*Sobs+j]=0.0;
      }
    }
  }
}

// Uniform on the log scale between -8 and 8
int allowable_rates(double *pc){
  int i, allowable=1;
  for (i=0;i<csnr;i++) {
    if ((pc[i]<exp(-8)) || (pc[i]>exp(8))) {
      allowable=0;
    }
  } 

  return(allowable);
}
int allowable_obsvar(int Sobs, double *pobsvar) {
  int i, allowable=1;
  for (i=0;i<Sobs;i++) {
    if ((pobsvar[i*Sobs+i]<-exp(20)) || (pobsvar[i*Sobs]>exp(3))) {
      allowable=0;
    }
  }
  return(allowable);
}

double log_prior_rates(double *pc) {
  return 0;
}

double log_prior_obsvar(int Sobs, double *pobsvar) {
  double lp=0;
  int i;

  for (i=0;i<Sobs;i++) {
    lp -= log(1+pobsvar[i*Sobs+i]*pobsvar[i*Sobs+i]); // Cauchy
  }

  return(lp); 
}

void get_jumpvar(int Sobs, double *pjv) {
  int i,j;

  // Zero first
  for (i=0;i<csnr+Sobs;i++) {
    for (j=0;j<csnr+Sobs;j++) {
	pjv[i*(Sobs+csnr)+j]=0.0;
    }
  }
  printf("%d %d\n",csnr,Sobs);
  for (i=0;i<csnr;i++) {
    pjv[i*(csnr+Sobs)+i]=0.01;
  }
}



#endif


void free_fixed_parms(void *pparms) {
  if (pparms != 0) {
    free(pparms);
  }
}

void get_SR(int *pS, int *pR) {
  *pS=csns; 
  *pR=csnr; 
}
void get_nits_and_seed(int *pnits, int *pseed) {
  *pnits=nits;
  *pseed=seed;
}

char *pget_system() {
  return psystem;
}
char *pget_obsfname() {
  return pobsfname;
}
char *pget_outfname() {
  return poutfname;
}

void get_hybrid_info(double *pN, double *phstar, double *pepsilon, 
		     double *pdeltat, int *pdo_sk){
  *pN=hybrid_N;
  *phstar=hybrid_hstar;
  *pepsilon=hybrid_epsilon;
  *pdeltat=hybrid_deltat;
  *pdo_sk=hybrid_use_sk;
}
