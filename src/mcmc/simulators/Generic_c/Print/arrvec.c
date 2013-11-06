#include <stdio.h>

void printdvec(int n,double *pv) {
  int i;
  for (i=0;i<n;i++) {
    printf("%f ",pv[i]);
  }
  printf("\n");
}
void fprintdvec(int n,double *pv, FILE *pfp) {
  int i;
  for (i=0;i<n;i++) {
    fprintf(pfp, "%f ",pv[i]);
  }
  fprintf(pfp,"\n");
}
void fprintdvecnoret(int n,double *pv, FILE *pfp) {
  int i;
  for (i=0;i<n;i++) {
    fprintf(pfp, "%f ",pv[i]);
  }
}

void printdmat(int nrow, int ncol,double *pv) {
  int i,j;
  for (i=0;i<nrow;i++) {
    for (j=0;j<ncol;j++) {
      printf("%f ",pv[i*ncol+j]);
    }
    printf("\n");
  }
  printf("\n");
}
void fprintdmat(int nrow, int ncol, double *pv, FILE *pfp) {
  int i,j;
  for (i=0;i<nrow;i++) {
    for (j=0;j<ncol;j++) {
      fprintf(pfp, "%f ",pv[i*ncol+j]);
    }
    fprintf(pfp,"\n");
  }
  fprintf(pfp,"\n");
}

void printdmat_precise(int nrow, int ncol,double *pv) {
  int i,j;
  for (i=0;i<nrow;i++) {
    for (j=0;j<ncol;j++) {
      printf("%15.14e ",pv[i*ncol+j]);
    }
    printf("\n");
  }
  printf("\n");
}

void printdarray(int nblock, int nrow, int ncol,double *pv) {
  int i,j,k;
  for (i=0;i<nblock;i++) {
    for (j=0;j<nrow;j++) {
      for (k=0;k<ncol;k++) {
	printf("%f ",pv[i*ncol*nrow+j*ncol+k]);
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");
}

/*********************************************/

void printivec(int n,int *pv) {
  int i;
  if (n>0) {
    for (i=0;i<n;i++) {
      printf("%i ",pv[i]);
    }
  }
  printf("\n");
}
void fprintivec(int n, int *pv, FILE *pfp) {
  int i;
  for (i=0;i<n;i++) {
    fprintf(pfp, "%i ",pv[i]);
  }
  fprintf(pfp,"\n");
}

void printimat(int nrow, int ncol,int *pv) {
  int i,j;
  if ((nrow>0) && (ncol>0)) {
    for (i=0;i<nrow;i++) {
      for (j=0;j<ncol;j++) {
	printf("%i ",pv[i*ncol+j]);
      }
      printf("\n");
    }
  }
  printf("\n");
}
void fprintimat(int nrow, int ncol, int *pv, FILE *pfp) {
  int i,j;
  for (i=0;i<nrow;i++) {
    for (j=0;j<ncol;j++) {
      fprintf(pfp, "%i ",pv[i*ncol+j]);
    }
    printf("\n");
  }
  fprintf(pfp,"\n");
}

void printstmat(int nrow, int ncol,size_t *pv) {
  int i,j;
  if ((nrow>0) && (ncol>0)) {
    for (i=0;i<nrow;i++) {
      for (j=0;j<ncol;j++) {
	printf("%i ",(int)pv[i*ncol+j]);
      }
      printf("\n");
    }
  }
  printf("\n");
}


void printiarr(int nblock, int nrow, int ncol,int *pv) {
  int i,j,k;
  if (nblock*nrow*ncol>0) {
    for (i=0;i<nblock;i++) {
      for (j=0;j<nrow;j++) {
	for (k=0;k<ncol;k++) {
	  printf("%i ",pv[i*ncol*nrow+j*ncol+k]);
	}
	printf("\n");
      }
      printf("\n");
    }
  }
  printf("\n");
}
