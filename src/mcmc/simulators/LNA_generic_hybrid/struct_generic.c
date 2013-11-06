#include <stdio.h>
#include <stdlib.h>

#include "struct_generic.h"
#include "hybrid.h"


void PrintReactionInfo(ReactionInfo *pI) {
  int S=pI->S;
  int R=pI->R;
  
  printf("\nREACTION INFO\n");
  printf("%d species in %d reactions.\n\n",S,R);
  printf("Rate constants: ");
  printdvec(R,pI->pc);

  printf("\nNet effects matrix:\n");
  printimat(R,S,pI->pA);

  printf("maxnreact %d; maxnspec %d\n",pI->maxnreact,pI->maxnspec);

  printf("\nSpecies affected (net) by reaction r:\n");
  printivec(R,pI->pNS);
  printimat(R,pI->maxnspec,pI->pS);

  printf("\nReactions that give a net change to species s:\n");
  printivec(S,pI->pNR);
  printimat(S,pI->maxnreact,pI->pR);

  printf("\nReactions that give a net change to species i and j:\n");
  printimat(S,S,pI->pNRR);
  printiarr(S,S,pI->maxnreact,pI->pRR);
}



void PrintSDEInfo(SDEInfo *pSI) {
  int S=pSI->S;
  int R=pSI->R;

  printf("\nSDEInfo\n%d species in %d reactions.\n",S,R);
  printf("\nReaction rates: ");
  printdvec(R,pSI->ph);
  printf("\nmu: ");
  printdvec(S,pSI->pmu);
  printf("\nS:\n");
  printdmat(S,S,pSI->pS);
  printf("\nS2:\n");
  printdmat(S,S,pSI->pS2);
  printf("\nF:\n");
  printdmat(S,S,pSI->pF);
  printf("\nHaV:\n");
  printdvec(S,pSI->pHaV);
  printf("\nHS2V:\n");
  printdmat(S,S,pSI->pHS2V);
  PrintHybridInfo(pSI->pHI);
}

void PrintLNARealisation(LNARealisation *pLR) {
  int S=pLR->S;
  printf("\nLNARealisation\n%d species at time %f.\n",S,pLR->t);
  printf("\n y: ");
  printdvec(S,pLR->py);
  printf("\n M: ");
  printdvec(S,pLR->pM); // discrepancy from deterministic true value 
  printf("\n x: ");
  printdvec(S,pLR->px); // true value 
}

void PrintObservation(Observation *pOI) {
  int S=pOI->S;
  int Sobs=pOI->Sobs;
  printf("\nObservation\n %d obs about %d species at time %f.\n",Sobs,S,pOI->t);
  printf("\n obs: ");
  printdvec(Sobs,pOI->pobs);
  //  printf("\n map P st obs = P X + error %p:\n",pOI->pP);
  //printdmat(Sobs,S,pOI->pP);
}
void PrintMeanVar(MeanVar *pMV) {
  int S=pMV->S;
  printf("\nMean and Variance at time %f, for %d species.\n",pMV->t,S);
  printf("\n mu: ");
  printdvec(S,pMV->pmu);
  printdmat(S,S,pMV->pSigma);
}

void PrintHybridInfo(HybridInfo *pHI) {
  int S=pHI->S;
  int R=pHI->R;
  
  printf("\nHYBRID INFO\n%d species in %d reactions (s=%d,f=%d).\n",S,R,pHI->nslow,pHI->nfast);
  printf("N=%f (S&K and GSG), epsilon=%f (GSG), hstar=%f (S&K), deltat (GSG) %f, SK? %d\n",pHI->N,pHI->epsilon,pHI->hstar,pHI->deltat,pHI->do_sk);
  printf("\nFast reation indicator: ");
  printivec(R,pHI->pifast);
  printf("\nSlow reaction indicator: ");
  printivec(R,pHI->pislow);
  printf("\nSlow rates: ");
  printdvec(R,pHI->pslow_rates);
  printf("\nSlow c: ");
  printdvec(R,pHI->pslow_c);
  printf("\nbstar: ");
  printdvec(S,pHI->pbstar);
  printf("\nbstarmax: ");
  printdvec(S,pHI->pbstarmax);
  printf("tot_slow=%f\n",pHI->tot_slow);
  printf("max_lamslow=%f\n",pHI->max_lamslow);
}


void PrintLNARealisationArray(int n, LNARealisation **ppRI) {
  int i;
  printf("\nLNARealisation Array\n\n");
  for (i=0;i<n;i++) {
    printf("\n**%d**\n",i);
    PrintLNARealisation(ppRI[i]);
  }
}

void PrintObservationArray(int n, Observation **ppOI) {
  int i;
  printf("\nObservation Array\n\n");
  for (i=0;i<n;i++) {
    printf("\n**%d**\n",i);
    PrintObservation(ppOI[i]);
  }
}

void PrintMeanVarArray(int n, MeanVar **ppMV) {
  int i;
  printf("\nMeanVar Array\n\n");
  for (i=0;i<n;i++) {
    printf("\n**%d**\n",i);
    PrintMeanVar(ppMV[i]);
  }
}


void FPrintRealisation(LNARealisation *pLR, FILE *pfp) {
  fprintf(pfp,"%f %d ",pLR->t,pLR->S);
  //  fprintdvecnoret(pLR->S,pLR->py,pfp);
  fprintdvec(pLR->S,pLR->px,pfp);
}

void FPrintAllRealisations(int n, LNARealisation **ppRI, FILE *pfp) {
  int i;
  fprintf(pfp,"%d\n",n);
  for (i=0;i<n;i++) {
    FPrintRealisation(ppRI[i],pfp);
  }
}

// ********** Allocate and Free ************


void AllocSDEInfo(ReactionInfo *pI, SDEInfo *pSI) {
  int S=pI->S;
  pSI->S=S;
  pSI->R=pI->R;
  pSI->ph=(double*) malloc(pSI->R*sizeof(double));
  pSI->pmu=(double*) malloc(S*sizeof(double));
  pSI->pF=(double*) malloc(S*S*sizeof(double));
  pSI->pS=(double*) malloc(S*S*sizeof(double));
  pSI->pS2=(double*) malloc(S*S*sizeof(double));
  pSI->pHaV=(double*) malloc(S*sizeof(double));
  pSI->pHS2V=(double*) malloc(S*S*sizeof(double));
  pSI->pscratch=(double*) malloc((S*7+S*S*20+pSI->R)*sizeof(double));
  //  printf("Allocated %d words to pSI->pscratch\n",S*7+S*S*20+pSI->R);
  pSI->pRI=pI;
  pSI->pHI=(HybridInfo*) malloc(sizeof(HybridInfo));
  AllocHybridInfo(pI,pSI->pHI);
}

void DeAllocSDEInfo(SDEInfo *pSI) {
  free(pSI->ph);
  free(pSI->pmu);
  free(pSI->pF);
  free(pSI->pS);
  free(pSI->pS2);
  free(pSI->pHaV);
  free(pSI->pHS2V);
  free(pSI->pscratch);
  DeAllocHybridInfo(pSI->pHI);
  free(pSI->pHI);
}

void AllocLNARealisation(ReactionInfo *pRI, LNARealisation *pLR) {
  int S=pRI->S;
  pLR->S=S;
  //  pLR->t=pLI->t;
  
  pLR->py = (double*)malloc(S*sizeof(double));
  pLR->pM= (double*)malloc(S*sizeof(double)); //discrep from determin value 
  pLR->px= (double*)malloc(S*sizeof(double)); //true value 
}

void DeAllocLNARealisation(LNARealisation *pLR) {
  free(pLR->py);
  free(pLR->px);
  free(pLR->pM);
}

void AllocObservation(int S, int Sobs, Observation *pOI, double *pP) {
  pOI->S=S;
  pOI->Sobs=Sobs;
  pOI->pobs = (double*)malloc(Sobs*sizeof(double));
  pOI->pP=pP;
}

void DeAllocObservation(Observation *pOI) {
  free(pOI->pobs);
}


void AllocMeanVar(ReactionInfo *pRI, MeanVar *pMV) {
  int S=pRI->S;
  pMV->S=S;
  pMV->pmu = (double*)malloc(S*sizeof(double));
  pMV->pSigma = (double*)malloc(S*S*sizeof(double));
}

void DeAllocMeanVar(MeanVar *pMV) {
  free(pMV->pmu);
  free(pMV->pSigma);
}


void AllocLNARealisationArray(int n_times, ReactionInfo *pRI, 
			      LNARealisation **ppLR) {
  int i;
  for (i=0;i<n_times; i++) {
    ppLR[i]=(LNARealisation*)malloc(sizeof(LNARealisation));
    AllocLNARealisation(pRI,ppLR[i]);
  }
}
void AllocObservationArray(int n_times, int S, int Sobs, Observation **ppOI,
			   double *pP) {
  int i;
  for (i=0;i<n_times; i++) {
    ppOI[i]=(Observation*)malloc(sizeof(Observation));
    AllocObservation(S,Sobs,ppOI[i],pP);
  }
}
void AllocMeanVarArray(int n_times, ReactionInfo *pRI, MeanVar **ppMV) {
  int i;
  for (i=0;i<n_times; i++) {
    ppMV[i]=(MeanVar*)malloc(sizeof(MeanVar));
    AllocMeanVar(pRI,ppMV[i]);
  }
}

void DeAllocLNARealisationArray(int n_times, LNARealisation **ppLR) {
  int i;
  for (i=0;i<n_times; i++) {
    DeAllocLNARealisation(ppLR[i]);
    free(ppLR[i]);
  }
}

void DeAllocObservationArray(int n_times, Observation **ppOI) {
  int i;
  for (i=0;i<n_times; i++) {
    DeAllocObservation(ppOI[i]);
    free(ppOI[i]);
  }
}
void DeAllocMeanVarArray(int n_times, MeanVar **ppMV) {
  int i;
  for (i=0;i<n_times; i++) {
    DeAllocMeanVar(ppMV[i]);
    free(ppMV[i]);
  }
}


void AllocHybridInfo(ReactionInfo *pI, HybridInfo *pHI) {
  pHI->S=pI->S;
  pHI->R=pI->R;
  pHI->pifast=(int*) malloc(pHI->R*sizeof(int));
  pHI->pislow=(int*) malloc(pHI->R*sizeof(int));
  //  pHI->pfast_rates=(double*) malloc(pHI->R*sizeof(double));
  pHI->pslow_rates=(double*) malloc(pHI->R*sizeof(double));
  pHI->pslow_c=(double*) malloc(pHI->R*sizeof(double));
  pHI->pbstar=(double*) malloc(pHI->S*sizeof(double));
  pHI->pbstarmax=(double*) malloc(pHI->S*sizeof(double));
  get_hybrid_info(&pHI->N,&pHI->hstar,&pHI->epsilon,&pHI->deltat,&pHI->do_sk);
}

void DeAllocHybridInfo(HybridInfo *pHI) {
  free(pHI->pislow);
  free(pHI->pifast);
  //  free(pHI->pfast_rates);
  free(pHI->pslow_rates);
  free(pHI->pslow_c);
  free(pHI->pbstar);
  free(pHI->pbstarmax);
}
