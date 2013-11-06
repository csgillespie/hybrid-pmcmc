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

//#define autoreg_1234_k
//#define LotkaVolterra
//#define GG
#define GGfriendly
//#define CSTEST
//#define CSTEST2

static void ForceNonNeg(double *pv, int n) {
  int i;
  for (i=0;i<n;i++) {
    pv[i]*=(pv[i]>0);
  }
}

static int hybrid_use_sk=0; // 0= use new decision method, 1=S&K
static double hybrid_hstar=4; // S&K

// COLIN - change this one - it is N*/\epsilon* of equation (10)
// I think N* of 10 and epsilon* of 0.25 is the minimal reasonable
// this gives hybrid_N=40
static double hybrid_N=60; // CS and S&K

// COLIN - change this to alter the frequency of checks on which reactions are 
// fast and which are slow
static double hybrid_deltat=0.1; 

// Condition on epsilon and Deltat is irrelevant for LNA as deltat stiff
static double hybrid_epsilon=0.25; // CS - irrelevant

#ifdef autoreg_1234_k

// 5 species and 8 reactions
// 0=DNA, 1=r=RNA, 2=P, 3=P2, 4=DNA.P2

// but   DNA.P2=k-DNA with k known
// 0 DNA + P2 -> DNA.P2
// 1 DNA.P2 -> DNA + P2
// 2 DNA -> DNA + RNA
// 3 RNA-> RNA + P
// 4 2P->P2
// 5 P2->2P
// 6 RNA->0
// 7 P->0

static int csns=4;
static int csnr=8;
static int nits=2000;
static int seed=1232123; 
//static int seed=1234; 
static double k=10;
static char psystem[]="autoreg 4 species and k\n";

//char pobsfname[]="../../Obs/GW2005Data/GW50x4.obs";
//char pobsfname[]="../../Obs/GW2005Data/GW100x4.obs";
//char pobsfname[]="../../Obs/GW2005Data/GW500x4.obs";

// My datasets

static char pobsfname[]="../../Obs/auto50x4.obs100";
//static char pobsfname[]="../../Obs/auto50x3.obs100";
//static char pobsfname[]="../../Obs/auto50x4ge.obs100";
//static char pobsfname[]="../../Obs/auto50x3ge.obs100";
//static char pobsfname[]="../../Obs/auto50x4pe.obs100";

static char poutfname[]="Output/auto50x4.mcmc100";
//static char poutfname[]="Output/auto50x3.mcmc100";
//static char poutfname[]="Output/auto50x4ge.mcmc100";
//static char poutfname[]="Output/auto50x3ge.mcmc100";
//static char poutfname[]="Output/auto50x4pe.mcmc100";
//static char poutfname[]="Output/nothing.mcmc";


void get_A(int *pA, int S, int R) {
  int pA0[]={-1,0,0,-1};
  int pA1[]={1,0,0,1};
  int pA2[]={0,1,0,0};
  int pA3[]={0,0,1,0};
  int pA4[]={0,0,-2,1};
  int pA5[]={0,0,2,-1};
  int pA6[]={0,-1,0,0};
  int pA7[]={0,0,-1,0};

  int i;

  assert(S==csns);
  assert(R==csnr);

  for (i=0;i<csns;i++) {
    pA[0*csns+i]=pA0[i];
    pA[1*csns+i]=pA1[i];
    pA[2*csns+i]=pA2[i];
    pA[3*csns+i]=pA3[i];
    pA[4*csns+i]=pA4[i];
    pA[5*csns+i]=pA5[i];
    pA[6*csns+i]=pA6[i];
    pA[7*csns+i]=pA7[i];
  }
}

void constrain(double *px) {
  if (px[0]<0) { // too little DNA so let reaction 1 happen
    px[3] -= px[0]; 
    px[0] =0;
  }
  else if (px[0]>k) { // too much DNA so let reaction 0 happen
    px[3] -= (px[0]-k);
    px[0] = k;
  }
  if (px[3]<0) { // too little P2 so let reaction 4 happen
    px[2] += 2*px[3];
    px[3]=0;
  }
  if (px[1]<0) { // too little RNA so let reaction 2 happen or reverse 6
    px[1] = 0;
  }
  if (px[2]<1) { // too little P so let reaction 3 happen or reverse 7
    px[2]=1;
  }
}

void get_h(double *ph, double *px, double *pc, void *pparms) {
  ph[0]=pc[0]*px[0]*px[3];
  ph[1]=pc[1]*(k-px[0]);
  ph[2]=pc[2]*px[0];
  ph[3]=pc[3]*px[1];
  ph[4]=pc[4]*0.5*px[2]*(px[2]-1);
  ph[5]=pc[5]*px[3];
  ph[6]=pc[6]*px[1];
  ph[7]=pc[7]*px[2];
}

void get_dhdx(double *pdhdx, double *px, double *pc, void *pparms) {
  vec_set(pdhdx,csns*csnr,0.0);

  pdhdx[0*csns+0]=pc[0]*px[3];
  pdhdx[0*csns+3]=pc[0]*px[0];

  pdhdx[1*csns+0]=-pc[1];


  pdhdx[2*csns+0]=pc[2];

  pdhdx[3*csns+1]=pc[3];

  pdhdx[4*csns+2]=pc[4]*(px[2]-0.5);
  
  pdhdx[5*csns+3]=pc[5];

  pdhdx[6*csns+1]=pc[6];

  pdhdx[7*csns+2]=pc[7];
}

void get_d2hdxmdxn(double *pd2hdxmdxn, double *px, double *pc, void *pparms,int hcpt) {
  vec_set(pd2hdxmdxn,16,0.0);

  switch(hcpt) {
  case 0: {
    pd2hdxmdxn[0*4+3]=pc[0];
    pd2hdxmdxn[3*4+0]=pc[0];
    break;
  }
  case 4: {
    pd2hdxmdxn[2*4+2]=pc[4];
    break;
  }

  }
}
void* get_fixed_parms() {
  double *pk=(double*)malloc(sizeof(double));
  *pk=k;
  return(pk);
}


void get_prior_for_X0(double *pmean, double *pvar) {
  int i,j;

  for (i=0;i<csns;i++) {
    pmean[i]=5; 
    for (j=0;j<csns;j++) {
      if (i==j) {
	pvar[i*csns+j]=100;
      }
      else {
	pvar[i*csns+j]=0.0;
      }
    }
  }
}

void initialise_rates(double *pc) {
  int i;
  // G&W(2005) values
  pc[0]=0.1; pc[1]=0.7; pc[2]=0.35; pc[3]=0.2;
  pc[4]=0.1; pc[5]=0.9; pc[6]=.3; pc[7]=.1;
  //pc[0]=10;pc[1]=70; // CSTMP 
  //pc[0]=-2.3; pc[1]=-.26; pc[2]=-.62; pc[3]=-5.31;
  //pc[4]=-3.3; pc[5]=-.92; pc[6]=.45; pc[7]=-8;
  //for (i=0;i<8;i++) {
  //  pc[i]=exp(pc[i]);
  // }


}

void initialise_obsvar(int Sobs, double *pobsvar) {
  int i,j;
  for (i=0;i<Sobs;i++) {
    for (j=0;j<Sobs;j++) {
      if (i==j) {
	pobsvar[i*Sobs+j]=1;
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
  double l=exp(-10), u=exp(8);
  for (i=0;i<csnr;i++) {
    if ((pc[i]<l) || (pc[i]>u)) {
      allowable=0;
    }
  } 

  return(allowable);
}
int allowable_obsvar(int Sobs, double *pobsvar) {
  int i, allowable=1;
  double l=exp(-20), u=exp(3);
  for (i=0;i<Sobs;i++) {
    if ((pobsvar[i*Sobs+i]<l) || (pobsvar[i*Sobs]>u)) {
      allowable=0;
    }
  }
  return(allowable);
}

double log_prior_rates(double *pc) {
  double lp=0;
  int i;
  for (i=0;i<csnr;i++) {
    lp-= log(1+4*pc[i]*pc[i]); // Cauchy with mode (on log scale) at c=1/2
    lp +=log(pc[i]);
  }
  return lp;
}

double log_prior_obsvar(int Sobs, double *pobsvar) {
  double lp=0;
  int i;
  
  for (i=0;i<Sobs;i++) {
    double v= pobsvar[i*Sobs+i];
    lp -= log(1+v*v); // Cauchy
    lp += log(v);
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
  
  //#define simple
#ifdef simple
  for (i=0;i<csnr+Sobs;i++) {
    pjv[i*(Sobs+csnr)+i]=.02;
  }
#endif
  
//#define auto_A // estimates based on 50x4.obs100
#ifdef auto_A
  pjv[0]=0.9;pjv[1]=0.8;
  pjv[13]=0.9;
  pjv[26]=0.05; pjv[30]=0.04;
  pjv[39]=0.05; pjv[43]=0.02;
  pjv[52]=.9;pjv[53]=0.8;
  pjv[65]=.9;
  pjv[78]=0.07;
  pjv[91]=0.05;
#endif

  //#define auto_B // estimates based on 50x3.obs100
#ifdef auto_B
  pjv[0]=1.6;pjv[1]=0.7;pjv[2]=0.5;
  pjv[12]=0.9;
  pjv[24]=0.4; 
  pjv[36]=0.5;pjv[40]=0.4
  pjv[48]=.9;pjv[49]=0.8;
  pjv[60]=.9;
  pjv[72]=0.07;
  pjv[84]=0.4;
#endif

  //#define auto_C
#ifdef auto_C  // estimates based on auto50x4ge.obs100
  pjv[0]=1.0;pjv[1]=0.9;
  pjv[13]=1.0;
  pjv[26]=0.4; pjv[30]=0.4; pjv[33]=-0.14;
  pjv[39]=0.4; pjv[43]=0.3;
  pjv[52]=.8;pjv[53]=0.7;
  pjv[65]=.8;
  pjv[78]=0.5;pjv[81]=-0.14;
  pjv[91]=0.3;
  pjv[104]=0.26;
  pjv[117]=0.15;
  pjv[130]=1.0;
  pjv[143]=0.14;
#endif

  //#define auto_D // estimates based on 50x3ge.obs100
#ifdef auto_D
  pjv[0]=1.8;pjv[1]=0.7;pjv[2]=0.7;
  pjv[12]=1.0;//pjv[17]=0.34;
  pjv[24]=0.8;//pjv[32]=0.15;
  pjv[36]=0.5;pjv[40]=0.4;
  pjv[48]=1.1;pjv[49]=1.0;
  pjv[60]=1.1;
  pjv[72]=0.4;
  pjv[84]=0.36;
  pjv[96]=0.08;
  pjv[108]=0.7;
  pjv[120]=0.21;
#endif

#define auto_E
#ifdef auto_E  // estimates based on auto50x4pe.obs100
  pjv[0]=1.2;pjv[1]=1.1;
  pjv[13]=1.2;
  pjv[26]=0.37; pjv[30]=0.35; pjv[33]=-0.14;
  pjv[39]=0.30; pjv[43]=0.29;// pjv[47]=-0.10;
  pjv[52]=0.9;pjv[53]=0.8;
  pjv[65]=0.9;
  pjv[78]=0.37;pjv[81]=-0.14;
  pjv[91]=0.30;//pjv[95]=-0.10;
  pjv[104]=0.07;
  pjv[117]=0.17;
  pjv[130]=0.10;
  pjv[143]=0.10;
#endif


#ifndef simple
  for (i=1;i<csnr;i++) {
    for (j=0;j<i;j++) {
      pjv[i*(Sobs+csnr)+j]=pjv[j*(Sobs+csnr)+i];
    }
  }
  vec_el_op(pjv,pjv,(Sobs+csnr)*(Sobs+csnr),2.38/(double)csnr/2.5,'*');
  //  vec_el_op(pjv,pjv,(Sobs+csnr)*(Sobs+csnr),2.38/csnr/2.0,'*');
  //vec_el_op(pjv,pjv,(Sobs+csnr)*(Sobs+csnr),2.38/csnr/5.0,'*');
#endif  

  printf("Jump variance matrix:\n");
  printdmat(csnr+Sobs,csnr+Sobs,pjv);
}



#endif

#ifdef LotkaVolterra

// 2 species and 3 reactions
// 0=Pred, 1=Prey

// but   DNA.P2=k-DNA with k known
// 0 Pred + Prey -> 2 x Pred
// 1 Pred ->  Pred -1
// 2 Prey -> Prey + 1

static int csns=2;
static int csnr=3;
static int nits=105000;
static int seed=1232123; 
//static int seed=1234; 
static char psystem[]="Lotka-Volterra\n";

//static char pobsfname[]="../../Obs/LV/LV50x2ne.obs101";
static char pobsfname[]="../../Obs/LV/LV50x1ne.obs101";

//static char poutfname[]="Output/LVtest.mcmc";

//static char poutfname[]="Output/LV50x2ne.mcmc101";
static char poutfname[]="Output/LV50x1ne.mcmc101";
//static char poutfname[]="Output/LV50x2ne.mcmc101K";
//static char poutfname[]="Output/LV50x1ne.mcmc101K";


void get_A(int *pA, int S, int R) {
  int pA0[]={1,-1};
  int pA1[]={-1,0};
  int pA2[]={0,1};

  int i;

  assert(S==csns);
  assert(R==csnr);

  for (i=0;i<csns;i++) {
    pA[0*csns+i]=pA0[i];
    pA[1*csns+i]=pA1[i];
    pA[2*csns+i]=pA2[i];
  }
}

void constrain(double *px) {
  if (px[0]<0) { 
    px[0] =0;
  }
  if (px[1]<0) { 
    px[1]=0;
  }
}

void get_h(double *ph, double *px, double *pc, void *pparms) {
  ph[0]=pc[0]*px[0]*px[1];
  ph[1]=pc[1]*px[0];
  ph[2]=pc[2]*px[1];
}

void get_dhdx(double *pdhdx, double *px, double *pc, void *pparms) {
  vec_set(pdhdx,csns*csnr,0.0);

  pdhdx[0*csns+0]=pc[0]*px[1];
  pdhdx[0*csns+1]=pc[0]*px[0];

  pdhdx[1*csns+0]=pc[1];

  pdhdx[2*csns+1]=pc[2];
}

void get_d2hdxmdxn(double *pd2hdxmdxn, double *px, double *pc, void *pparms,int hcpt) {
  vec_set(pd2hdxmdxn,16,0.0);

  switch(hcpt) {
  case 0: {
    pd2hdxmdxn[0*csns+1]=pc[0];
    pd2hdxmdxn[1*csns+0]=pc[0];
    break;
  }
  }
}
void* get_fixed_parms() {
  return(0);
}

void get_prior_for_X0(double *pmean, double *pvar) {
  int i,j;

  pmean[0]=20;
  pmean[1]=80;

  pvar[0]=0.000001;
  pvar[1]=0;
  pvar[2]=0;
  pvar[3]=0.000001;
}

void initialise_rates(double *pc) {
  int i;
  pc[0]=0.01; pc[1]=0.6; pc[2]=0.3; 
}

void initialise_obsvar(int Sobs, double *pobsvar) {
  int i,j;
  for (i=0;i<Sobs;i++) {
    for (j=0;j<Sobs;j++) {
      if (i==j) {
	pobsvar[i*Sobs+j]=0.00001;
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
  double l=exp(-10), u=exp(8);
  for (i=0;i<csnr;i++) {
    if ((pc[i]<l) || (pc[i]>u)) {
      allowable=0;
    }
  } 

  return(allowable);
}
int allowable_obsvar(int Sobs, double *pobsvar) {
  int i, allowable=1;
  double l=exp(-20), u=exp(3);
  for (i=0;i<Sobs;i++) {
    if ((pobsvar[i*Sobs+i]<l) || (pobsvar[i*Sobs]>u)) {
      allowable=0;
    }
  }
  return(allowable);
}

double log_prior_rates(double *pc) {
  double lp=0;
  int i;
  for (i=0;i<csnr;i++) {
    lp-= log(1+10*sqrt(pc[i]*pc[i])); 
    lp +=0.5*log(pc[i]);
  }
  return lp;
}

double log_prior_obsvar(int Sobs, double *pobsvar) {
  double lp=0;
  int i;
  
  for (i=0;i<Sobs;i++) {
    double v= pobsvar[i*Sobs+i];
    lp -= log(1+v*v); // Cauchy
    lp += log(v);
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
  
  for (i=0;i<csnr;i++) {
    pjv[i*(Sobs+csnr)+i]=.02; // LNAKF .02  KOM .012
  }
  
  printf("Jump variance matrix:\n");
  printdmat(csnr+Sobs,csnr+Sobs,pjv);
}



#endif // LV


#ifdef GG
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


#ifdef GGfriendly // only difference is that final reaction produces 21 X2
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





#ifdef CSTEST
// 1 species and 2 reactions
static int csns=1;
static int csnr=2;
static char psystem[] = "CS Test System\n";
static char pobsfname[]="../../Obs/cstest4x1.obs";

void get_A(int *pA, int S, int R) {
  int pA0[]={1};
  int pA1[]={-1};

  int i;

  assert(S==csns);
  assert(R==csnr);

  for (i=0;i<csns;i++) {
    pA[0*csns+i]=pA0[i];
    pA[1*csns+i]=pA1[i];
  }
}

void get_h(double *ph, double *px, double *pc, void *pparms) {
  ph[0]=pc[0];
  ph[1]=pc[1]*px[0];
  ForceNonNeg(ph,csnr);
}

void get_dhdx(double *pdhdx, double *px, double *pc, void *pparms) {
  pdhdx[0*csns+0]=0; 
  pdhdx[1*csns+0]=pc[1];
}

void* get_fixed_parms() {
  return 0;
}

void PriorForX0(MeanVar *pprior) {
  int i,j;
  int S=pprior->S;
  pprior->t=0.0;

  for (i=0;i<S;i++) {
    pprior->pmu[i]=100;
    for (j=0;j<S;j++) {
      if (i==j) {
	pprior->pSigma[i*S+j]=2500;
      }
      else {
	pprior->pSigma[i*S+j]=0.0;
      }
    }
  }
}

void initialise_rate_constants(double *ptheta) {
  double *pc=ptheta;
  int i;
  // So that equilibrium is at 100 molecules
  pc[0]=10; pc[1]=0.1; 
  for (i=0;i<csnr;i++) {
    pc[i]=log(pc[i]);
  }
}

void initialise_observation_variances(double *ptheta,int Sobs) {
  double *pvar=ptheta+csnr;
  int i;
  for (i=0;i<Sobs;i++) {
    pvar[i]=log(1);
  }
}



#endif


#ifdef CSTEST2
// 3 species and 6 reactions
static int csns=3;
static int csnr=6;
char pfobsfname[]="../../Obs/cstest4x3.obs";
char pfobsfname[]="../../Obs/cstest4x2.obs";
char pfobsfname[]="../../Obs/cstest500x3.obs";

void get_A(int *pA, int S, int R) {
  int pA0[]={1,0,0};
  int pA1[]={-1,0,0};
  int pA2[]={0,1,0};
  int pA3[]={0,-1,0};
  int pA4[]={0,0,1};
  int pA5[]={0,0,-1};

  int i;

  assert(S==csns);
  assert(R==csnr);

  for (i=0;i<csns;i++) {
    pA[0*csns+i]=pA0[i];
    pA[1*csns+i]=pA1[i];
    pA[2*csns+i]=pA2[i];
    pA[3*csns+i]=pA3[i];
    pA[4*csns+i]=pA4[i];
    pA[5*csns+i]=pA5[i];
  }
}

void get_h(double *ph, double *px, double *pc, void *pparms) {
  ph[0]=pc[0];
  ph[1]=pc[1]*px[0];
  ph[2]=pc[2];
  ph[3]=pc[3]*px[1];
  ph[4]=pc[4];
  ph[5]=pc[5]*px[2];
  ForceNonNeg(ph,csnr);
}

void get_dhdx(double *pdhdx, double *px, double *pc, void *pparms) {
  int i,j;
  for (i=0;i<csnr;i++) {
    for (j=0;j<csns;j++) {
      pdhdx[i*csns+j]=0.0;
    }
  }
  pdhdx[1*csns+0]=pc[1]; 
  pdhdx[3*csns+1]=pc[3]; 
  pdhdx[5*csns+2]=pc[5];
}

void* get_fixed_parms() {
  return 0;
}

void PriorForX0(MeanVar *pprior) {
  int i,j;
  int S=pprior->S;
  pprior->t=0.0;

  for (i=0;i<S;i++) {
    pprior->pmu[i]=150;
    for (j=0;j<S;j++) {
      if (i==j) {
	pprior->pSigma[i*S+j]=400;
      }
      else {
	pprior->pSigma[i*S+j]=0.0;
      }
    }
  }
}

void initialise_rate_constants(double *ptheta) {
  double *pc=ptheta;
  int i;
  // So that equilibrium is at 100 molecules
  pc[0]=10; pc[1]=0.1; pc[2]=15; pc[3]=0.15; pc[4]=5; pc[5]=0.05; 
  for (i=0;i<csnr;i++) {
    pc[i]=log(pc[i]);
  }
}

void initialise_observation_variances(double *ptheta,int Sobs) {
  double *pvar=ptheta+csnr;
  int i;
  for (i=0;i<Sobs;i++) {
    pvar[i]=log(0.01);
  }
}


void get_jumpvar(double *pjv, int Sobs) {
  int i,j;

  for (i=0;i<csnr+Sobs;i++) {
    for (j=0;j<csnr+Sobs;j++) {
      if ((i==j) && (i<csnr)) {
	pjv[i*(Sobs+csnr)+j]=0.003/5; // CSTMP just added for TEST
      }
      else {
	pjv[i*(Sobs+csnr)+j]=0.0;
      }
    }
  }

  //  for (i=csnr;i<csnr+Sobs;i++) {
  // pjv[i*(Sobs+csnr)+i]=0.01;
  //}

  printf("Jump variance matrix:\n");
  printdmat(csnr+Sobs,csnr+Sobs,pjv);
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
