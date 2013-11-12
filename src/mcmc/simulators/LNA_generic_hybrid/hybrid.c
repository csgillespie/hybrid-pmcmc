#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "struct_generic.h"
#include "hybrid.h"

static int nfast2=0;
static int nslow2=0;
static int nfast4=0;
static int nslow4=0;
static int nr2=0;
static double totx2r2=0;
static double tothr2=0;
static double tothdtr2=0;
static double totepsilon=0;
//static double totS=0;


void PrintNFastSlow() {
  double dr2=(double)nr2;
  printf("%d fast2 and %d slow2\n",nfast2,nslow2);
  printf("%d fast4 and %d slow4\n",nfast4,nslow4);
  nfast2=0;nslow2=0;
  nfast4=0;nslow4=0;
  printf("mean h %f, mean hdt %f, mean x2 %f, mean eps %f\n",tothr2/dr2,tothdtr2/dr2,totx2r2/dr2, totepsilon/dr2);
  totx2r2=0; tothr2=0; tothdtr2=0; totepsilon=0; nr2=0;
}

#ifdef SKPOSS
static int FastReaction(double *px,int *pA, double hdt, int S, 
			double N, double hstar, double epsilon, int use_sk) {
  int fast=1, i=0;

  if (hdt<1) {
    hdt=1; // GSG
  }

  if (use_sk) {
    fast=fast* (hdt > hstar); // S&K
  }
  while (fast && (i<S)) {
    fast=fast*(int)(px[i]>= N*fabs((double)pA[i])); // S&K and GSG
    // Commented out code below unnecessary as stiff solver takes care of it 
    if (!use_sk) {
      fast=fast*(int)(px[i]>= hdt*fabs((double)pA[i])/epsilon); // GSG
     }
    i++;
  }
  // fast=0; // CSTMP all slow for testing
  // fast=1; // CSTMP all fast for testing
  return(fast);
}
#endif

static int FastReaction(double *px,int *pA, double hdt, int S, 
			double N, double hstar, double epsilon, int use_sk) {
  int fast=1, i=0;
  double a;

  if (hdt<1) {
    hdt=1; // GSG
  }
 
  while (fast && (i<S)) {
    a=fabs((double)pA[i]);
    fast=fast*(int)(px[i]>= N*a); // S&K and GSG
    fast=fast*(int)(px[i]>= hdt*a/epsilon); // GSG
    i++;
  }
  // fast=0; // CSTMP all slow for testing
  // fast=1; // CSTMP all fast for testing
  return(fast);
}

void SetHybridInfo(HybridInfo *pHI, int *pA, double *px, double *ph, 
		   double *pc,double deltat) {
  int i, R=pHI->R, S=pHI->S;
  pHI->nslow=0;
  pHI->nfast=0;

  for (i=0;i<R;i++) {
    pHI->pifast[i]=FastReaction(px,pA+i*S,ph[i]*deltat,S,
				pHI->N,pHI->hstar,pHI->epsilon,pHI->do_sk);
    pHI->pislow[i]=1-pHI->pifast[i];
    pHI->nfast+=pHI->pifast[i];
    pHI->nslow+=pHI->pislow[i];
    pHI->pslow_c[i] = pc[i] * (double) pHI->pislow[i];
  }
  for (i=0;i<S;i++) {
    pHI->pbstar[i]=0;
    pHI->pbstarmax[i]=0;
  }

  pHI->max_lamslow=0; // can always zero at start of a new interval
  pHI->tot_slow=0;
}

void SetHybridRates(HybridInfo *pHI, double *ph) {
  int i;
  int R=pHI->R;
  pHI->tot_slow=0.0;

  for (i=0;i<R;i++) {
    pHI->pslow_rates[i]=(double) pHI->pislow[i] * ph[i];
    pHI->tot_slow += pHI->pslow_rates[i];
  }
  if (pHI->tot_slow>pHI->max_lamslow) {
    pHI->max_lamslow = pHI->tot_slow; // running maximum
  }
}
void HybridModifyRates(double *ph, HybridInfo *pHI) {
  int i, R=pHI->R;

  for (i=0;i<R;i++) {
    ph[i]=(double) pHI->pifast[i] * ph[i];
  }
}

void HybridModifyDerivs(double *pdhdx, HybridInfo *pHI) {
  int i,j;
  int R=pHI->R, S=pHI->S;

  for (i=0;i<R;i++) {
    if (pHI->pislow[i]) {
      for (j=0;j<S;j++) {
	pdhdx[i*S+j]=0.0;
      }
    }
  }
}

int *pISlow(HybridInfo *pHI) {
  return(pHI->pislow);
}

double *pHybridSlowRates(HybridInfo *pHI) {
  return(pHI->pslow_rates);
}


int nSlow(HybridInfo *pHI) {
  return(pHI->nslow);
}
double TotSlowRate(HybridInfo *pHI) {
  return(pHI->tot_slow);
}

int IChooseSlowReaction(HybridInfo *pHI, double uniform) {
  int i=0,ichoose=-1;
  double run_tot_slow=0;
  double mx=uniform*pHI->tot_slow;
  while (run_tot_slow<mx) {
    if (pHI->pislow[i]) {
      ichoose=i;
      run_tot_slow += pHI->pslow_rates[i];
    }
    i++;
  }
  
  return(ichoose); // -1 if no slow reactions + avoid double arithmetic probs
}

