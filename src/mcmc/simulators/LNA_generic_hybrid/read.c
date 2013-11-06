#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "struct_generic.h"
#include "../Generic_c/Print/arrvec.h"


void ReadEffectsMatrices(char *pfname, int **ppL, int **ppR, int *pS, int *pR) {
  FILE *pfp=fopen(pfname,"r");
  char *pline = (char*) malloc(301); // max line length
  char *p;
  int i,j,S,R;
  int currpos=0;

  p=fgets(pline,300,pfp);
  while (pline[currpos]==' ') {
    currpos++;
  }
  *pS=(int)(pline[currpos++]-'0');
  while (pline[currpos]==' ') {
    currpos++;
  }
  *pR=(int)(pline[currpos]-'0');

  *ppL= (int *) malloc(*pS * (*pR) * sizeof(int));
  *ppR= (int *) malloc(*pS * (*pR) * sizeof(int));

  printf("S and R are %d %d\n",*pS,*pR);

  S=*pS; R=*pR;
  for (i=0;i<R;i++) {
    currpos=0;

    p=fgets(pline,300,pfp);
    printf("%s\n",pline);

    for (j=0;j<S;j++) {
      while (pline[currpos]==' ') {
	currpos++;
      }
      (*ppL)[i*S+j]=(int)(pline[currpos++]-'0');
    }
    for (j=0;j<S;j++) {
      while (pline[currpos]==' ') {
	currpos++;
      }
      (*ppR)[i*S+j]=(int)(pline[currpos++]-'0');
    }
  }

  fclose(pfp);
}

// Read Observation file

void ReadNObsAndSobs(FILE *pfp,int *pnobs, int *pobslen) {
  fscanf(pfp,"%d %d",pnobs,pobslen);
}

void BracketNext(char *pline, int *pi, int *pj) {
  *pi=*pj+1;
  while (((pline[*pi]<'0') || (pline[*pi]>'9')) && (pline[*pi]!='.')) {
    (*pi)++;
  }
  *pj=*pi;
  while (((pline[*pj]>='0') && (pline[*pj]<='9')) || (pline[*pj]=='.')) {
    (*pj)++;
  }
}

// Fill the Sobs x S matrix that obs = P X + error
void ReadP(FILE *pfp, int maxlen, int S, int Sobs, double *pP) {
  char *p, *pline = (char*) malloc(maxlen+1);
  int i,j,m,n;

  p = fgets(pline, maxlen, pfp); // ignore remains of the first line 
  p = fgets(pline, maxlen, pfp); // ignore blank line

  for (m=0;m<Sobs;m++) {
    p = fgets(pline, maxlen, pfp);
    j=-1;
    for (n=0;n<S;n++) {
      BracketNext(pline,&i,&j);
      pline[j]=0;
      pP[m*S+n] = atof(pline+i);
    }
  }
  free(pline);
}

void LineToObs(char *pline, Observation *pOI) {
  int j=-1;
  int i,k;

  BracketNext(pline,&i,&j);
  pline[j]=0;
  pOI->t = atof(pline+i);

  // Obs files used to have Sobs in each line - removed.
  //  BracketNext(pline,&i,&j);
  //pline[j]=0;
  //assert(pOI->Sobs==atoi(pline+i));

  for (k=0;k<pOI->Sobs;k++) {
    BracketNext(pline,&i,&j);
    pline[j]=0;
    pOI->pobs[k]=atof(pline+i);
  }
}

void ReadAllObs(FILE *pfp, int maxlen, int nobs, Observation **ppOI) {
  char *p, *pline = (char*) malloc(maxlen+1);
  int i,j;

  p = fgets(pline, maxlen, pfp); //  blank line

  for (i=0;i<nobs;i++) {
    p = fgets(pline, maxlen, pfp);
    LineToObs(pline,ppOI[i]);
  }
  free(pline);
}

//int main(int argc, char *argv[]) {
    void irrelevant(int x) {
  Observation **ppOI;
  char pinfile[]="autoregfull.obs";
  FILE *pfp=fopen(pinfile,"r");
  double *pP;
  int nobs;
  int S=5,Sobs;
  int maxlen=200;

  if (pfp == NULL) {
    printf("Observation file <%s> does not exist.\n",pinfile);
    exit(1);
  }
  
  ReadNObsAndSobs(pfp,&nobs,&Sobs);
  printf("nobs %d, Sobs %d\n",nobs,Sobs);
  ppOI=(Observation **)malloc(nobs*sizeof(Observation *));
  pP=(double*)malloc(S*Sobs*sizeof(double));
  ReadP(pfp,maxlen,S,Sobs,pP);

  AllocObservationArray(nobs, S, Sobs, ppOI, pP);

  ReadAllObs(pfp, maxlen, nobs, ppOI);

  PrintObservationArray(nobs,ppOI);
  DeAllocObservationArray(nobs,ppOI);
  free(ppOI);
  free(pP);
  fclose(pfp);
}
