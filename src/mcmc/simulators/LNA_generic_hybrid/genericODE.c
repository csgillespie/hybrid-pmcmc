#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <time.h>
#include <sys/timeb.h>

#include "../Generic_c/LinearAlgebra/cslinalg.h"
#include "../Generic_c/Print/arrvec.h"
#include "struct_generic.h"
#include "hybrid.h"
#include "read.h"
#include "system_spec.h"
#include "genericODE.h"
#include "random.h"


// Which equations should the ODE solver solve?
// It always performs a deterministic integral for the mean

//#define ODE_VARIANCE    // Simpler and faster than generator
//#define TWO_MOMENT    // Add 2MA approximation terms to the mean and VARIANCE 
#define ODE_GENERATOR // Use the generator method
#define ODE_TAU       // Required for hybrid method - needs ODE_GENERATOR

// Allowable combinations
// ODE_VARIANCE                (for LNAKF)
// ODE_VARIANCE and TWO_MOMENT (for 2MAKF)
// ODE_GENERATOR               (for LNAKF - using generator)
// ODE_GENERATOR and ODE_TAU   (for hybrid systems)
 
//#define PRINT_derivinfo  // for debugging
//#define PRINT_integrate
//#define PRINT_HYBRID

//#define COMPARE_VARIANCES // if using both ODE_VARIANCE and ODE_GENERATOR

static double epsilon=0.001;

typedef struct _ode_info {
  int S;
  double t;
  int len_metay;
  int len_metay_with_sq;

  double *pmetay; // has pu, pG, pPsihalf, pVhalf, ptau, + pPsi, pV
  double *pu; // log of current deterministic value
  double *pV; // variance matrix
  double *pV_half;
  double *pG; // inverse of generator for ODE
  double *pPsi; // for covariance
  double *pPsi_half;
  double *ptau; // (i,j) = jth time for ith BM
  double *pm; // mean of the perturbation

  double *pGinv; // generator
} ODEInfo;

static double *puf(double *pmetay, int S) {
  return(pmetay);
}
static double *pGf(double *pmetay, int S) {
#ifdef ODE_GENERATOR
  return(pmetay+S);
#else
  return(0);
#endif
}
static double *pPsi_halff(double *pmetay, int S) {
#ifdef ODE_GENERATOR
  return(pmetay+S+S*S);
#else
  return(0);
#endif
}
static double *pPsif(double *pmetay, int S) {
#ifdef ODE_GENERATOR
  int index=S+S*S+S*(S+1)/2;
#ifdef ODE_VARIANCE
  index += S*(S+1)/2
#endif
#ifdef ODE_TAU
    index += S;
#endif
  return(pmetay+index);
#else
  return(0);
#endif
}
static double *pV_halff(double *pmetay, int S) {
#ifdef ODE_VARIANCE
  int index=S;
#ifdef ODE_GENERATOR
  index += S*S + S*(S+1)/2;
#endif
  return(pmetay+index);
#else
  return(0);
#endif
}
static double *pVf(double *pmetay, int S) {
#ifdef ODE_VARIANCE
  int index;
  index=S+S*(S+1)/2;
#ifdef ODE_GENERATOR
  index += 2*S*S + S*(S+1)/2;
#endif
#ifdef ODE_TAU
  index += S;
#endif
  return(pmetay+index);
#else
  return(0);
#endif
}
static double *ptauf(double *pmetay, int S) {
#ifdef ODE_TAU
  int index;
  index=S;
#ifdef ODE_GENERATOR
  index += S*S + S*(S+1)/2;
#else
  assert(5==7);
#endif
#ifdef ODE_VARIANCE
  index += S*(S+1)/2;
#endif
  return(pmetay+index);
#else
  return(0);
#endif
}
// ODE size
static int len_metay(int S) {
  int len = S;
#ifdef ODE_GENERATOR
  len += S*S + S*(S+1)/2;
#ifdef ODE_TAU
  len +=S;
#endif
#endif
#ifdef ODE_VARIANCE
  len += S*(S+1)/2;
#endif
  return(len);
}
// memory to alloc
static int len_metay_sq(int S) {
  int len = S;
#ifdef ODE_GENERATOR
  len += S*S + S*(S+1)/2 + S*S;
#ifdef ODE_TAU
  len +=S;
#endif
#endif
#ifdef ODE_VARIANCE
  len += S*(S+1)/2 + S*S;
#endif
  return(len);
}


// Convert between lower triangular matrix (n(n+1)/2 vals) and 
// square symmetric nxn matrix
static void LowerTriToSquare(double *pL, double *pS, int n) {
  int i,j,ind=0;

  for (i=0;i<n;i++) {
    for (j=0;j<=i;j++) {
      pS[i*n+j]=pL[ind];
      pS[j*n+i]=pL[ind];
      ind++;
    }
  }
  assert(ind==n*(n+1)/2);
}
static void SquareToLowerTri(double *pL, double *pS, int n) {
  int i,j,ind=0;

  for (i=0;i<n;i++) {
    for (j=0;j<=i;j++) {
      pL[ind]=pS[i*n+j];
      ind++;
    }
  }
  assert(ind==n*(n+1)/2);
}

static void yS_to_yL(double *pmetay, double *pmetay_sq, int S) {
#ifdef ODE_GENERATOR
  SquareToLowerTri(pPsi_halff(pmetay,S),pPsif(pmetay_sq,S),S);
#endif
#ifdef ODE_VARIANCE
  SquareToLowerTri(pV_halff(pmetay,S),pVf(pmetay_sq,S),S);
#endif
}

static void yL_to_yS(double *pmetay, double *pmetay_sq,int S) {
#ifdef ODE_GENERATOR
  LowerTriToSquare(pPsi_halff(pmetay,S),pPsif(pmetay_sq,S),S);
#endif
#ifdef ODE_VARIANCE
  LowerTriToSquare(pV_halff(pmetay,S),pVf(pmetay_sq,S),S);
#endif
}


void* pODEInfoAlloc(void *pRIv) {
  ReactionInfo *pI = (ReactionInfo*)pRIv;
  ODEInfo *pLI=(ODEInfo*)malloc(sizeof (ODEInfo));
  int S=pI->S;
  pLI->S=S;
  int index;
  //  pLI->t=time;
  
  pLI->len_metay = len_metay(S);
  pLI->len_metay_with_sq = len_metay_sq(S);

  // All variables for ODE solver in one contigous space
  pLI->pmetay = (double*)malloc(pLI->len_metay_with_sq*sizeof(double));

  // But for ease of access...
  pLI->pu= pLI->pmetay; // current deterministic value
#ifdef ODE_GENERATOR
  pLI->pG= pGf(pLI->pmetay,S);
  pLI->pPsi_half= pPsi_halff(pLI->pmetay,S);
  pLI->pPsi= pPsif(pLI->pmetay,S);
#endif
#ifdef ODE_VARIANCE
  pLI->pVhalf=pVhalf(pLI->pmetay,S);
  pLI->pV=pVf(pLI->pmetay,S);
#endif
#ifdef ODE_TAU
  pLI->ptau= ptauf(pLI->pmetay,S);//vec of marginal BM times
#endif

  /* printf("Allocating %d words (ODEsize %d) for metay at %p,\npu %p, pG %p, pPsi_hal, %p, ptau %p, pVhalf %p, pPsi %p, pV %p\n",pLI->len_metay_with_sq,pLI->len_metay,pLI->pmetay,pLI->pu, pLI->pG,pLI->pPsi_half, pLI->ptau, pLI->pV_half,pLI->pPsi,pLI->pV); */

  // Also 
#ifdef ODE_GENERATOR
  pLI->pGinv= (double*)malloc(S*S*sizeof(double)); 
#endif
  return((void*) pLI);
}

void ODEInfoDeAlloc(void *pOI) {
  ODEInfo *pLI = (ODEInfo*) pOI;
  free(pLI->pmetay);
#ifdef ODE_GENERATOR
  free(pLI->pGinv);
#endif
  free(pLI);
}

void ODEInfoPrint(void *pOI) {
  ODEInfo *pLI = (ODEInfo*) pOI;
  int S=pLI->S;
  printf("\nODEInfo\n%d species at time %f.\n",S,pLI->t);
  printf("\n u: ");
  printdvec(S,pLI->pu);
#ifdef ODE_GENERATOR
  printf("\n G\n");
  printdmat(S,S,pLI->pG);
  printf("\n Ginv\n");
  printdmat(S,S,pLI->pGinv);
  printf("\n Psi\n");
  printdmat(S,S,pLI->pPsi);
  printf("\n Psi_half\n");
  printdvec(S*(S+1)/2,pLI->pPsi_half);
#endif
#ifdef ODE_VARIANCE
  printf("\n V\n");
  printdmat(S,S,pLI->pV);
  printf("\n V_half\n");
  printdvec(S*(S+1)/2,pLI->pV_half);
#endif
#ifdef ODE_TAU
  printf("\n tau\n");
  printdvec(S,pLI->ptau);
#endif
}

void ODEInfoArrayPrint(int n, void **ppOI) {
  ODEInfo **ppLI = (ODEInfo**) ppOI;
  int i;
  printf("\nODEInfo Array\n\n");
  for (i=0;i<n;i++) {
    printf("**%d**\n",i);
    ODEInfoPrint(ppLI[i]);
  }
}

void** ppODEInfoArrayAlloc(int n_times, ReactionInfo *pRI) {
  ODEInfo **ppLI = (ODEInfo**) malloc(n_times *sizeof(ODEInfo*));
  int i;
  for (i=0;i<n_times; i++) {
    ppLI[i]=(ODEInfo*)pODEInfoAlloc(pRI);
  }
  return((void**)ppLI);
}
void ODEInfoArrayDeAlloc(int n_times, void **ppOI) {
  ODEInfo **ppLI = (ODEInfo**) ppOI;
  int i;
  for (i=0;i<n_times; i++) {
    ODEInfoDeAlloc(ppLI[i]);
  }
  free(ppOI);
}

 
/*********************************************/


// [Number of] species that experience a net change in each reaction
static void NumberOfSpecies(int R, int S, int *pA, int *pNS) {
  int i,j;
  for (i=0;i<R;i++) {
    pNS[i]=0;
    for (j=0;j<S;j++) {
      if (pA[i*S+j]!=0) {
	(pNS[i]) += 1;
      }
    }
  }
}

static void Species(int R, int S, int *pA, int maxnspec, int *pS) {
  int i,j;
  for (i=0;i<R;i++) {
    for (j=0;j<maxnspec;j++) {
      pS[i*maxnspec+j]=-1;
    }
  }
  for (i=0;i<R;i++) {
    int s=0;
    for (j=0;j<S;j++) {
      if (pA[i*S+j]!=0) {
	pS[i*maxnspec+s] = j;
	s++;
      }
    }
  }
}

// [Number of] species that contribute to reaction rate (2P->P2 means two)
static void NumberOfLeftSpecies(int R, int S, int *pA, int *pNLS) {
  int i,j;
  for (i=0;i<R;i++) {
    pNLS[i]=0;
    for (j=0;j<S;j++) {
      if (pA[i*S+j]!=0) { 
	(pNLS[i]) -= pA[i*S+j];
      }
    }
    assert(pNLS[i]<=2);
  }
}

static void LeftSpecies(int R, int S, int *pA, int *pLS) {
  int i,j;
  for (i=0;i<R;i++) {
    for (j=0;j<2;j++) {
      pLS[i*2+j]=-1;
    }
  }
  for (i=0;i<R;i++) {
    int s=0;
    for (j=0;j<S;j++) {
      if (pA[i*S+j]!=0) {
	pLS[i*2+s] = j;
	s++;
      }
      if (pA[i*S+j]==2) {
	pLS[i*2+s] = j;
	s++;
      } 
    }
  }
}

// [Number of] reactions that cause a net change to each species
static void NumberOfReactionsNetChange(int R, int S, int *pA, int *pNR) {
  int i,j,tot;
  for (j=0;j<S;j++) {
    tot=0;
    for (i=0;i<R;i++) {
      if (pA[i*S+j] != 0) {
	tot++;
      }
    }
    pNR[j]=tot;
  }
}

static void ReactionsNetChange(int R, int S, int *pA, int maxnreact, int *pR) {
  int i,j;
  for (j=0;j<S;j++) {
    for (i=0;i<maxnreact;i++) {
      pR[j*maxnreact+i]=-1;
    }
  }

  for (j=0;j<S;j++) {
    int r=0;
    for (i=0;i<R;i++) {
      if (pA[i*S+j] != 0) {
	pR[j*maxnreact+r] =i;
	r++;
      }
      assert(r<=maxnreact);
    }
  }
}


static void NumberOfReactionsNetNetCross(int R, int S, int *pA, int *pNRR) {
  int i,j,k;
  for (i=0;i<S;i++) {
    for (j=0;j<S;j++) { 
      // Reactions involving species i and j
      pNRR[i*S+j]=0;
      for (k=0;k<R;k++) {
	if ((pA[k*S+i]!=0) && (pA[k*S+j]!=0)) {
	  pNRR[i*S+j] ++;
	}
      } 
    }
  }
}

static void ReactionsNetNetCross(int R, int S, int *pA, int maxnreact, int *pRR) {
  int i,j,k;
  for (i=0;i<S;i++) {
    for (j=0;j<S;j++) { 
      // Reactions involving species i and j
      for (k=0;k<maxnreact;k++) {
	pRR[i*S*maxnreact+j*maxnreact+k]=-1;
      }
    }
  }

  for (i=0;i<S;i++) {
    for (j=0;j<S;j++) { 
      // Reactions involving species i and j
      int r=0;
      for (k=0;k<R;k++) {
	if ((pA[k*S+i]!=0) && (pA[k*S+j]!=0)) {
	  pRR[i*S*maxnreact+j*maxnreact+r] =k;
	  r++;
	}
      } 
    }
  }
}


static void NumberOfReactionsNetLeftCross(int R, int S, int *pA, int *pAL, 
					  int *pNRRL) {
  int i,j,k;
  for (i=0;i<S;i++) {
    for (j=0;j<S;j++) { 
      // Reactions involving species i and j
      pNRRL[i*S+j]=0;
      for (k=0;k<R;k++) {
	if ((pA[k*S+i]!=0) && (pAL[k*S+j]!=0)) {
	  pNRRL[i*S+j] ++;
	}
      } 
    }
  }
}

static void ReactionsNetLeftCross(int R, int S, int *pA, int *pAL, 
				  int maxnreact, int *pRRL) {
  int i,j,k;

  for (i=0;i<S;i++) {
    for (j=0;j<S;j++) { 
      // Reactions involving species i and j
      for (k=0;k<maxnreact;k++) {
	pRRL[i*S*maxnreact+j*maxnreact+k]=-1;
      }
    }
  }


  for (i=0;i<S;i++) {
    for (j=0;j<S;j++) { 
      // Reactions involving species i and j
      int r=0;
      for (k=0;k<R;k++) {
	if ((pA[k*S+i]!=0) && (pAL[k*S+j]!=0)) {
	  pRRL[i*S*maxnreact+j*maxnreact+r] =k;
	  r++;
	}
      } 
    }
  }
}

// R reactions, S species, A = net effect matrix
void AllocReactionInfo(ReactionInfo *pRI) {
  int S,R,i,j;
  int *pA;

  get_SR(&pRI->S,&pRI->R);
  S=pRI->S; R=pRI->R;

  pRI->pA=(int*) malloc(R*S*sizeof(int)); // net effects
  pRI->pc = (double*) malloc(R*sizeof(double));
  pRI->pNS=(int*) malloc(R*sizeof(int)); //  number of species net changed by re
  pRI->pNR=(int*) malloc(S*sizeof(int)); // # reactions for each species
  pRI->pNRR=(int*) malloc(S*S*sizeof(int)); //#reactions for each pair of spec
  pRI->pparms=get_fixed_parms(); // mallocs and sets

  pA=pRI->pA;

  // Net Effects matrix and any additional parameters

  get_A(pA,S,R);

  // For fast slow decision need species affected net by each reaction

  NumberOfSpecies(R,S,pA,pRI->pNS);
  pRI->maxnspec=MaxOfIntVector(R,pRI->pNS);
  pRI->pS=(int*) malloc(R*(pRI->maxnspec)*sizeof(int)); //reactions for each sp
  Species(R,S,pA,pRI->maxnspec,pRI->pS);

  // Overall mean rate of change for species i depends on net effects matrix

  NumberOfReactionsNetChange(R,S,pA,pRI->pNR);
  pRI->maxnreact=MaxOfIntVector(S,pRI->pNR);
  pRI->pR=(int*) malloc(S*(pRI->maxnreact)*sizeof(int)); //reactions for each sp
  ReactionsNetChange(R,S,pA,pRI->maxnreact,pRI->pR);

  // For variance calculations need reactions that have net effect on species 
  // i and j

  NumberOfReactionsNetNetCross(R,S,pA,pRI->pNRR);
  pRI->pRR=(int*) malloc(S*S*(pRI->maxnreact)*sizeof(int));
  ReactionsNetNetCross(R,S,pA,pRI->maxnreact,pRI->pRR);
}

void DeAllocReactionInfo(ReactionInfo *pI) {
  free(pI->pc);
  free(pI->pA);
  free(pI->pNS);
  free(pI->pS);
  free(pI->pNR);
  free(pI->pR);
  free(pI->pNRR);
  free(pI->pRR);
  free_fixed_parms(pI->pparms);
}

//*****************************************************

// Deterministic SDE term
// Hardcode for GGfriendly
static void get_a(double *pmu, double *px, double *pcf) {
  double tmp=pcf[4]*px[0]*px[1];
  pmu[0]=pcf[0]-pcf[2]*px[0]-tmp;
  pmu[1]=pcf[1]-pcf[3]*px[1]+20*tmp;
}

// Matrix of first derivatives of reaction rates (for LNA)
// Hardcode for GGfriendly
static void get_F(double *pF, double *px, double *pcf) {
  pF[0]=-pcf[2]-pcf[4]*px[1];
  pF[1]=-pcf[4]*px[0];
  pF[2]=20*pcf[4]*px[1];
  pF[3]=-pcf[3]+20*pcf[4]*px[0];
}

// Variance for SDE
// Hardcode for GGfriendly 
static void get_S2(double *pS2, double *px, double *pcf) {
  double tmp=pcf[4]*px[0]*px[1];
  pS2[0]=pcf[0]+pcf[2]*px[0]+tmp;
  pS2[1]=20*tmp;
  pS2[2]=pS2[1];
  pS2[3]=pcf[1]+pcf[3]*px[1]+400*tmp;
}


static void SDEInfoFrom_y(double *pV, SDEInfo *pSI,
			  double *peta, HybridInfo *pHI) {
  ReactionInfo *pRI=pSI->pRI;
  double pcf[5]; // rate constants but zero if slow
  int i;

  for (i=0;i<5;i++) {
    if (pHI->pislow[i]) {
      pcf[i]=0;
    }
    else {
      pcf[i]=pRI->pc[i];
    }
  }

  get_h(pSI->ph,peta,pRI->pc,pRI->pparms);
  SetHybridRates(pHI,pSI->ph);
  //  printf("h\n");
  //printdvec(pSI->S,pSI->ph);
  //  HybridModifyRates(pSI->ph,pHI);
  //printdvec(pSI->S,pSI->ph);
  get_a(pSI->pmu,peta,pcf); // aka "a" - the multiplier of dt
  get_F(pSI->pF,peta,pcf);
  get_S2(pSI->pS2,peta,pcf);

  

#ifdef ODE_TAU
  matrix_square_root(pSI->S, pSI->pS2, pSI->pS);
#endif

#ifdef PRINT_derivinfo
  printf("eta:\n");
  printdvec(pSI->S,peta);
#ifdef ODE_VARIANCE
  printf("V:\n");
  printdmat(pSI->S,pSI->S,pV);
#endif
  printf("c:\n");
  printdvec(pRI->R,pRI->pc);
  printf("h\n");
  printdvec(pSI->R,pSI->ph);
  printf("mu\n");
  printdvec(pSI->S,pSI->pmu);
  printf("F\n");
  printdmat(pSI->S,pSI->S,pSI->pF);
  printf("S2:\n");
  printdmat(pSI->S,pSI->S,pSI->pS2);
#ifdef ODE_TAU
  printf("S:\n");
  printdmat(pSI->S,pSI->S,pSI->pS);
#endif
#endif
}

void EnsurePositive(double *pvec, int n) {
  int i;
  double epsilon=0.01;
  for (i=0;i<n;i++) {
    if (pvec[i]<epsilon) {
      pvec[i]=epsilon;
    }
  }
}

// ********** ODE solver *************

// If only one of the species is involved in fast reactions

// y = [u0,G00,psi00,tau0] 
// params is [cf0,cf2,cs0,cs1,cs2,cs3,cs4,xfixed,bmx,lammxsl,hybrid] 
// (c0 and c2 potentially fast) or
// y = [u1,G11,psi11,tau1]
// params is [cf1,cf3,cs1,cs0,cs3,cs2,cs4,xfixed,bmx,lammxsl,hybrid] 
// (c1 and c3 potentially fast)
extern void lsoda_derivs_special(int *podedim, double *pt, double *py, double *pdydt, double *pparm) {
  double u = py[0];
  double eta = exp(u)-epsilon;
  double G = py[1];
  double G2=G*G;
  double detadt = pparm[0]-pparm[1]*eta;
  double dudt = detadt/exp(u);
  double dGdt=-pparm[1]*G;
  double S2=pparm[0]+pparm[1]*eta;
  double dPsidt=S2/G2;
  double *pcs=pparm+2;
  double *pxfixed=pcs+5; // other cpt is fixed
  double *phybrid=pxfixed+4;
  double tmp;

  pdydt[0]=dudt;
  pdydt[1]=dGdt;
  pdydt[2]=dPsidt;

  //printf("dudt=%f\n",dudt);
  //printf("dGdt=%f\n",dGdt);
  //printf("dPsi=%f\n",dPsidt);

  if (*phybrid==1) {
    double dtaudt=S2/G2;
    double pbstar[2];
    double *pbmx=pxfixed+1;
    double *plammaxsl=pbmx+2;

    pdydt[3]=dtaudt;

  // Running maximum for bstar

    pbstar[0]=G*(pcs[2]+pcs[4]* *pxfixed); 
    pbstar[1]=pcs[3]+pcs[4]*eta;

    // printf("dtau=%f\n",dtaudt);
    //printf("b*=%f,%f\n",pbstar[0],pbstar[1]);

    tmp=fabs(pbstar[0]);
    if (tmp>pbmx[0]) {
      pbmx[0]=tmp;
    }
    tmp=fabs(pbstar[1]);
    if (tmp>pbmx[1]) {
      pbmx[1]=tmp;
    }

    //    printf("b*max=%f,%f\n",pbmx[0],pbmx[1]);

  // Running maximum for lam tot slow
    tmp=pcs[0]+pcs[1]+pcs[2]*eta+pcs[3]* *pxfixed+pcs[4]*eta* *pxfixed;
    if (tmp>*plammaxsl) {
      *plammaxsl=tmp;
    }

    //  printf("lammaxsl=%f\n",*plammaxsl);
  }
}

extern void lsoda_derivs(int *podedim, double *pt, double *py, double *pdydt, double *pparams) {
  SDEInfo *pSI= (SDEInfo*) pparams;
  HybridInfo *pHI=pSI->pHI;
  const int S=pSI->S;

  double *peta = pSI->pscratch;
  double *pdydt_sq = peta+S;
  double *ptmp = pdydt_sq+3*S+3*S*S+S*(S+1); // u+G+psi+V+tau+psih+Vh
  double *ptmp2 = ptmp+2*S+3*S*S;
  double *ptmp3= ptmp2 +S*S;
  double *pmetay=(double*)py;
  double *pu=puf(pmetay,S);
  double *pdudt=puf(pdydt,S);
#ifdef ODE_GENERATOR
  double *pG=pGf(pmetay,S);
  double *pPsi=pPsif(pmetay,S);
  double *pdGdt=pGf(pdydt,S);
  double *pdPsidt=pPsif(pdydt_sq,S);
#endif
  double *pV=pVf(pmetay,S);
#ifdef ODE_VARIANCE
  double *pdVdt=pVf(pdydt_sq,S);
#endif
#ifdef ODE_TAU
  double *ptau=ptauf(pmetay,S);
  double *pdtaudt=ptauf(pdydt,S);
  double *pS=pSI->pS;
#endif
  double *pF=pSI->pF;
  double *pS2=pSI->pS2;
  int i,j,k,pos;  
  
  pos=S;

  yL_to_yS((double*)py,(double*)py,S);

#ifdef COMPARE_VARIANCES
  matrix_invert(S,pG,ptmp);
  mat_mult(ptmp,pPsi,ptmp2,S,S,S,S,0,0);
  mat_mult(ptmp2,ptmp,ptmp3,S,S,S,S,0,1);
  printf("Variance comparison\n");
  printdmat(S,S,ptmp3);
  printdmat(S,S,pV);
#endif

  vec_el_op(pu,peta,S,0,'e'); 
  vec_el_op(peta,peta,S,epsilon,'-'); // eta = exp(y) - epsilon
  constrain(peta);

  SDEInfoFrom_y(pV,pSI,peta,pHI);

  // d eta / dt =  mu   so  d log(eta+e) / dt = mu / (eta +e)=mu/exp(y) 
  for (i=0;i<S;i++) {
    pdudt[i]=pSI->pmu[i]/exp(pu[i]); // 
  }

#ifdef ODE_GENERATOR
  // dG/dt = -GF  nb G is the inverse generator here
  mat_mult((double*)pG,(double*)pF,pdGdt,S,S,S,S,0,0);
  for (i=0;i<S*S;i++) {
    pdGdt[i] = -pdGdt[i];
  }

  // dPsi_{ij}/dt = [G S2 G^t]_{ij}
  mat_mult((double*)pG,(double*)pS2,ptmp,S,S,S,S,0,0);
  mat_mult(ptmp,(double*)pG,pdPsidt,S,S,S,S,0,1);
#endif

  //  if (0) {
  //  printf("dudt=%f,%f\n",pdudt[0],pdudt[1]);
  // printf("dGdt=%f,%f,%f,%f\n",pdGdt[0],pdGdt[1],pdGdt[2],pdGdt[3]);
  //printf("dPsi=%f,%f,%f,%f\n",pdPsidt[0],pdPsidt[1],pdPsidt[2],pdPsidt[3]);
  //printf("dtau=%f,%f\n",pdtaudt[0],pdtaudt[1]);
  //}

#ifdef ODE_VARIANCE
  // dV/dt = FV + VF^t +S2
  mat_mult((double*)pF,(double*)pV,ptmp,S,S,S,S,0,0);
  mat_add_transpose(S,ptmp,ptmp);
  vecvec_el_op((double*)pS2,ptmp,pdVdt,S*S,'+');
#endif

#ifdef ODE_TAU
  if (pHI->hybridrun==1) {
    // dtau_{i}/dt= sum_j [GS]_{ij}^2
    // So that tau_i is the marginal time for the ith BM 
    for (i=0;i<S;i++) {
      double sum=0;
      for (j=0;j<S;j++) {
	double el=0; // i,j element of GS
	for (k=0;k<S;k++) {
	  el += pG[i*S+k]*pS[k*S+j];
	}
	sum += el*el;
      }
      pdtaudt[i]= sum;
    }

    matrix_invert(S,pG,ptmp2);
    get_b(ptmp,peta,pHI->pslow_c);
    mat_mult(ptmp2,ptmp,pHI->pbstar,S,S,S,1,1,0);

    for (i=0;i<S;i++) {
      if (fabs(pHI->pbstar[i])>pHI->pbstarmax[i]) {
	pHI->pbstarmax[i]=fabs(pHI->pbstar[i]);
      }
    }
#ifdef PRINT_HYBRID
    printf("b: ");
    printdvec(S,ptmp);
    printf("G:\n");
    printdmat(S,S,ptmp2);
    printf("bstar: ");
    printdvec(S,pHI->pbstar);
    printf("bstarmax: ");
    printdvec(S,pHI->pbstarmax);
    printf("totlamslowmax: %f\n",pHI->max_lamslow);
#endif
  }
  // else {
    //    for (i=0;i<S;i++) {
    //pdtaudt[i]=0.0;
    //}
  //}
#endif

  yS_to_yL(pdydt,pdydt_sq,S); // take non-blank bits of pdydt_sq i.e. V and psi
}

static void MeanToLog(double *pu,double *pMean, int S) {
  int i;

  for (i=0;i<S;i++) {
    pu[i]=log(pMean[i]+epsilon);
  }
}

static void LogToMean(double *pMean, double *pu, int S) {
  int i;

  for (i=0;i<S;i++) {
    pMean[i]=exp(pu[i])-epsilon;
  }
}

// Takes from Posterior Mean and Variance to initialise the LI
// Takes the previous mean of the perturbation for use with Komorowski
// Also initialises G to identity and tau to 0
void ODEInitialise(void *pOI, double t, double *pPostMean, double *pPostVar, double *peta) {
  ODEInfo *pLI = (ODEInfo*) pOI;
  int S=pLI->S;
  int i,j;

  pLI->t=t;

  MeanToLog(pLI->pu,pPostMean,S);

#ifdef ODE_GENERATOR
  memcpy(pLI->pPsi,pPostVar,S*S*sizeof(double));
#endif
#ifdef ODE_VARIANCE
  memcpy(pLI->pV,pPostVar,S*S*sizeof(double));
#endif

#ifdef ODE_GENERATOR  
  for (i=0;i<S;i++) {
    for (j=0;j<S;j++) {
      if (i==j) {
	pLI->pG[S*i+j]=1;
      }
      else {
	pLI->pG[S*i+j]=0;
      }
    }
#ifdef ODE_TAU 
    pLI->ptau[i]=0;
#endif
  }
#endif

}


void ODEExtract(void *pOI,double *pPriorMean, double *pPriorVar, double *peta, 
		double *pt, double *ptmp) {
  ODEInfo *pLI = (ODEInfo*) pOI;
  int i,S=pLI->S;
  double *pGen = pLI->pGinv;

  *pt=pLI->t;

  // mu_{prior} = y_{LNA} + G M = y_{LNA} since M=0

  LogToMean(pPriorMean,pLI->pu,S);

#ifdef ODE_GENERATOR
  // Sigma_{prior} = G Psi_{LNA} G^t - calculation checked and works!
  mat_mult(pGen,pLI->pPsi,ptmp,S,S,S,S,0,0);
  mat_mult(ptmp,pGen,pPriorVar,S,S,S,S,0,1);
#ifdef COMPARE_VARIANCES
  printf("Var from generator\n");
  printdmat(S,S,pPriorVar);
#endif
#endif

#ifdef ODE_VARIANCE
  memcpy(pPriorVar,pLI->pV,S*S*sizeof(double));
#ifdef COMPARE_VARIANCES
  printf("Var from simpler\n");
  printdmat(S,S,pPriorVar);
#endif
  // Multiplying large pPsi by small G can lead to minor asymmetries
  // which become major when the Kriging (posterior) variance is calculated
  // So squash them now
  SymmetriseSquareMatrix(S,pPriorVar);
#endif
}

static double qnorm(double a) { 
  return(gsl_cdf_ugaussian_Pinv(a));
}

// Must use pmetay - not pmetay_sq
static double get_max_intensity(HybridInfo *pHI, double *ptau) {
  int S=pHI->S, i;
  double u;
  double epsilon = 1e-6; // probability of exceeding threshold
  double tot;
  double sum_u_bstar=0;

  for (i=1;i<S;i++) {
    u = -sqrt(ptau[i])*qnorm(epsilon/(4*S));
    sum_u_bstar += u*pHI->pbstarmax[i];
  }
  tot = sum_u_bstar+pHI->max_lamslow;
#ifdef PRINT_HYBRID
  printf("sum ub* %f, qn %f,maxdetlamslow %f,\n combined: %f\n",sum_u_bstar,qnorm(epsilon/(4*S)),pHI->max_lamslow,tot);
#endif
  return(tot);
}

static void SampleAfterFast(ODEInfo *pLI, double *psample, double *pscratch, 
			    gsl_rng *r) {
  double tt;
  int S=pLI->S;
  double *pmean = pscratch;
  double *pvar = pscratch + S;
  double *ptmp = pvar + S*S;

  ODEExtract(pLI, pmean, pvar, ptmp, &tt, ptmp);
  SimulateTruncMVN(S,pmean,pvar,ptmp ,psample,r);
}

static void UpdateBySlow(HybridInfo *pHI, int *pA, double *pstate, 
			 gsl_rng *r, double *pscratch) {
  int rreact=IChooseSlowReaction(pHI,gsl_ran_flat(r,0.0,1.0));
  int i;

  for (i=0;i<pHI->S;i++) {
    pstate[i]+=(double) pA[rreact*pHI->S+i];
    if (pstate[i]<0) {
      pstate[i]=0;
    }
  }
#ifdef PRINT_HYBRID
	printf("Slow reaction was %d and new state vec is\n",rreact);
	printdvec(pHI->S,pstate);
#endif
}

// Integrate forward for time t - alters py, pG and ptau in ODEInfo
void ODEIntegrateSimulateOne(double t, SDEInfo *pSI, void *pOI, double *pzero,
			     double *pxcurr, gsl_rng *r, double *pt_int, 
			     double *pt_all_slow) {
  int i,S=pSI->S,R=pSI->R;
  ODEInfo *pLI = (ODEInfo*) pOI;
  HybridInfo *pHI=pSI->pHI;
  int odesize=pLI->len_metay;
  int odesize_notau=odesize-S;
  int shodesize=4;
  int shodesize_notau=3;
  double ttodo=t, tstart;
  double next_potential_hit=ttodo;
  double lammaxslow=-1, lamtotslow;
  double rerr=1.0e-4, aerr=1.0e-4;
  int actual_hit=1;
  double *pscratch= pSI->pscratch;
  double *ph=(double*) malloc(R*sizeof(double));
  double *pdhdx=(double*) malloc(R*S*sizeof(double));
  double *pxcurr_loop_start=(double*) malloc(S*sizeof(double));
  double integrate_time;

  *pt_int=0;
  *pt_all_slow=0;

  // Start with pxcurr for deterministic
  while (ttodo>0) {
    //    printf("ttdodo %f\n",ttodo);
    if (ttodo<pHI->deltat) {
      integrate_time=ttodo;
    }
    else {
      integrate_time=pHI->deltat;
    }

    // Hybrid: decide on fast/slow status for each reaction
    // only need to do this every actual hit or fixed hybrid time
    get_h(ph,pxcurr,pSI->pRI->pc,pSI->pRI->pparms);
    //    printf("pxcurr %f, %f\n",pxcurr[0],pxcurr[1]); // CSTMP
    SetHybridInfo(pHI,pSI->pRI->pA,pxcurr,ph,pSI->pRI->pc,integrate_time);
#ifdef PRINT_HYBRID
    printf("Starting loop: ttodo=%f, old next_potential_hit=%f, integrate_time=%f\n Fast/Slow decision (no rates yet):\n",ttodo,next_potential_hit,integrate_time);
    PrintHybridInfo(pHI);
    printdvec(S,pxcurr);
#endif
    if (pHI->nfast==0) { // All Slow
      SetHybridRates(pHI,ph);
      lamtotslow=TotSlowRate(pHI);
      next_potential_hit = -log(gsl_ran_flat(r,0.0,1.0))/lamtotslow;
      if (next_potential_hit < ttodo) { // not hybrid_time here as all slow
	*pt_all_slow += next_potential_hit;
	UpdateBySlow(pHI,pSI->pRI->pA,pxcurr,r,pscratch);
#ifdef PRINT_HYBRID
	printf("Slow reaction at %f as all slow (with integrate_time %f, ttodo %f)?\n",next_potential_hit,integrate_time, ttodo);
#endif
      }
      else {
	*pt_all_slow += ttodo;
      }
      ttodo -= next_potential_hit; // goes -ve when out of time

    }
    else { // At least one fast
      // Special hardcode for when at most two fast reactions and these
      // only change one species
      int shortcut=0;
      memcpy(pxcurr_loop_start,pxcurr,S*sizeof(double));

      if (pHI->pislow[4]) {
	if (pHI->pislow[1] && pHI->pislow[3]) {
	  shortcut=1;
	}
	else if (pHI->pislow[0] && pHI->pislow[2]) {
	  shortcut=2;
	}
      }
      //                        if (0) {
      if (shortcut>0) {
	double xfast, cfA, cfB;
	double pparm[11];
	double *pcs=pparm+2;
	double *pxfixed=pcs+5; // other cpt is fixed
	double *pbmx=pxfixed+1;
	double *plammaxsl=pbmx+2;
	double *phybrid=plammaxsl+1;
	double py[4];
	double *pc=pSI->pRI->pc;
	double ph[5],phs[5];
	double mn,sd,tots,tmp;
	double ptau[2];

	get_h(ph,pxcurr,pc,0);
	tots=0;
	for (i=0;i<5;i++) {
	  pcs[i]=pc[i]*(double)pHI->pislow[i];
	  tots+=ph[i]*(double)pHI->pislow[i];
	}
	*plammaxsl=tots; //  initialise lammaxslow
	pbmx[0]=0; pbmx[1]=0; // initialise bmx


	if (shortcut==1) {
	  xfast=pxcurr[0];
	  *pxfixed=pxcurr[1];
	  cfA=pc[0]*(1-pHI->pislow[0]); cfB=pc[2]*(1-pHI->pislow[2]);
	}
	else {
	  xfast=pxcurr[1];
	  *pxfixed=pxcurr[0];
	  cfA=pc[1]*(1-pHI->pislow[1]); cfB=pc[3]*(1-pHI->pislow[3]);
	  tmp=pcs[0]; pcs[0]=pcs[1]; pcs[1]=tmp;
	  tmp=pcs[2]; pcs[2]=pcs[3]; pcs[3]=tmp;
	}
	py[0]=log(xfast+epsilon);
	py[1]=1; // G
	py[2]=0; // psi
	py[3]=0; // tau
	pparm[0]=cfA; pparm[1]=cfB;

	*phybrid=1;
	tstart=0;
	VG_integrate(py, &tstart, &integrate_time, pparm, &shodesize, 
		     &lsoda_derivs_special, &rerr, &aerr);
	//	assert(7==8);

	if (shortcut==1) {
	  pHI->pbstarmax[0]=pbmx[0]; pHI->pbstarmax[1]=pbmx[1];
	  ptau[0]=py[3];ptau[1]=0;
	}
	else {
	  pHI->pbstarmax[1]=pbmx[0]; pHI->pbstarmax[0]=pbmx[1];
	  ptau[1]=py[3];ptau[0]=0;
	}
	pHI->max_lamslow=*plammaxsl;
	if (pHI->nslow>0) {
	  lammaxslow = get_max_intensity(pHI,ptau);
	  next_potential_hit = -log(gsl_ran_flat(r,0.0,1.0))/lammaxslow;

	  if (next_potential_hit < integrate_time) {
	    *phybrid=2;
	    tstart=0;
	    py[0]=log(xfast+epsilon);
	    py[1]=1; // G
	    py[2]=0; // psi
	    py[3]=0; // tau
	    (*pt_int)+=next_potential_hit;
	    VG_integrate(py, &tstart, &next_potential_hit, pparm, 
			 &shodesize_notau, 
			 &lsoda_derivs_special, &rerr, &aerr);
	  }
	}
	// Sim. the new fast x from 
	mn=exp(py[0])-epsilon;
	sd=py[1]*sqrt(py[2]);
	if (shortcut==1) {
	  pxcurr[0]=mn+sd*gsl_ran_gaussian(r,1.0);
	}
	else { // shortcut=2
	  pxcurr[1]=mn+sd*gsl_ran_gaussian(r,1.0);
	}
	
      }
      else {

	ODEInitialise(pOI,0,pxcurr,pzero,pzero); // set square;tau=0 and psi=0, G=I
	yS_to_yL(pLI->pmetay,pLI->pmetay,S); // square -> lower triangle
	get_dhdx(pdhdx,pxcurr,pSI->pRI->pc,pSI->pRI->pparms);

#ifdef PRINT_integrate
	printf("Before integration\n");
	ODEInfoPrint(pLI);
#endif
#ifdef PRINT_HYBRID
	printf("Before first integration of interval %f \n",ttodo);
	ODEInfoPrint(pLI);
#endif


	pHI->hybridrun=1; // collect maxima and taus
	tstart=0;
	(*pt_int)+=integrate_time;
	VG_integrate(pLI->pmetay, &tstart, &integrate_time, (double*)pSI, &odesize, 
		   &lsoda_derivs, &rerr, &aerr);
	//		assert(7==8);
#ifdef PRINT_HYBRID
	printf("Done first integration: 0 to %f\nNot made square.\n",integrate_time);
	ODEInfoPrint(pLI);
#endif

	if (pHI->nslow>0) {  // Mixture of fast and slow
	  lammaxslow = get_max_intensity(pHI,ptauf(pLI->pmetay,S));
	  next_potential_hit = -log(gsl_ran_flat(r,0.0,1.0))/lammaxslow;
#ifdef PRINT_HYBRID
	  printf("max intensity %f, new next_pot_hit %f, integrate_time %f\n",lammaxslow,next_potential_hit,integrate_time);
#endif

	  if (next_potential_hit < integrate_time) {
	    pHI->hybridrun=2;    // don't bother with maxima or taus
	    tstart=0;
	    ODEInitialise(pOI,0,pxcurr_loop_start,pzero,pzero); //tau=0, psi=0, G=I
	    yS_to_yL(pLI->pmetay,pLI->pmetay,S); // square -> lower triangle
	    (*pt_int)+=next_potential_hit;
	    VG_integrate(pLI->pmetay, &tstart, &next_potential_hit, (double*)pSI, &odesize_notau, &lsoda_derivs, &rerr, &aerr);
#ifdef PRINT_HYBRID
	    printf("Done second integration: 0 to %f\nNot made square.\n",next_potential_hit);
	    ODEInfoPrint(pLI);	
#endif
	  }
	} // end mixture of fast and slow
	yL_to_yS(pLI->pmetay,pLI->pmetay,S); // square is current now

#ifdef ODE_GENERATOR
	matrix_invert(S,pLI->pG,pLI->pGinv);
	SymmetriseSquareMatrix(S,pLI->pPsi); // remove minor asymmetries
#endif
      // pxcurr is current again from here

	SampleAfterFast(pLI, pxcurr, pscratch, r); // from truncated MVN; needs Ys
      } // end of choice with fast hardcode
      if ((pHI->nslow>0) && (next_potential_hit < integrate_time)) {
#ifdef PRINT_HYBRID
	printf("After 2nd integration (but before slow) state and ODE Info is \n");
	printdvec(S,pxcurr);
	ODEInfoPrint(pLI);
#endif
	get_h(ph,pxcurr,pSI->pRI->pc,pSI->pRI->pparms);

	SetHybridRates(pHI,ph);
	lamtotslow=TotSlowRate(pHI);
	if (gsl_ran_flat(r,0.0,1.0)*lammaxslow < lamtotslow) {
	  UpdateBySlow(pHI,pSI->pRI->pA,pxcurr,r,pscratch);
#ifdef PRINT_HYBRID
	  printf("Realised slow potential reaction at %f (with integrate_time %f)? %f/%f=%f.\n",next_potential_hit,integrate_time,lamtotslow,lammaxslow,lamtotslow/lammaxslow);
#endif
	}
	else {
#ifdef PRINT_HYBRID
      printf("Not realised slow potential reaction at %f (after integrate_time %f).\n",next_potential_hit,integrate_time);
      PrintHybridInfo(pHI);
#endif
	}
	ttodo -= next_potential_hit;// whether or not a hit	
      }
      else {
	ttodo -= integrate_time; 
      }

    } // end of at least one fast

  }     // end loop with pxcurr current

  pLI->t += t;

  //  PrintNFastSlow();

#ifdef PRINT_integrate
  printf("After integration\n");
  ODEInfoPrint(pLI);
#endif

  free(ph);
  free(pdhdx);
  free(pxcurr_loop_start);
}
