#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include "../Print/arrvec.h"

//******* Matrix routines ***********

void lu_invert(int n, double *pmx, double *pinv) {
  double *ptmpmx = (double*) malloc(n*n*sizeof(double));
  double *ptmpvc = (double*) malloc(n*sizeof(double));
  gsl_matrix_view m = gsl_matrix_view_array (ptmpmx, n, n);
  gsl_vector_view b  = gsl_vector_view_array (ptmpvc, n);
  gsl_permutation *ptmpp = gsl_permutation_alloc(n);
  gsl_vector *px = gsl_vector_alloc (n);
  int i,s,j;

  memcpy(ptmpmx,pmx,n*n*sizeof(double));
  gsl_linalg_LU_decomp(&m.matrix, ptmpp, &s);

  *ptmpvc=1;
  for (i=1;i<n;i++) {
    ptmpvc[i]=0;
  }
  gsl_linalg_LU_solve(&m.matrix, ptmpp, &b.vector, px);

  for (i=0;i<n;i++) {
    pinv[i*n]=px->data[i];
  }

  for (j=1;j<n;j++) {
    ptmpvc[j-1]=0;
    ptmpvc[j]=1;

    gsl_linalg_LU_solve(&m.matrix, ptmpp, &b.vector, px);

    for (i=0;i<n;i++) {
      pinv[i*n+j]=px->data[i];
    }
  }

  free(ptmpmx);
  free(ptmpvc);
  gsl_vector_free(px);
  gsl_permutation_free (ptmpp);
}

void matrix_invert(int n, double *pmx, double *pinv) {
  if (n==1) {
    *pinv=1/ *pmx;
  }
  else if (n==2) {
    double det=(pmx[0]*pmx[3]-pmx[2]*pmx[1]);
    pinv[0]=pmx[3]/det; pinv[1]=-pmx[1]/det; 
    pinv[2]=-pmx[2]/det; pinv[3]=pmx[0]/det;
  }
  else {
    lu_invert(n,pmx,pinv);
  }
}


// Cholesky decompose 2x2 matrix [a,b\\  b,c]
void chol22(double a, double b, double c,
	    double *pL11, double *pL21, double *pL22) {
  if (a == 0) {
    *pL11=0;
    *pL21=0;
    *pL22=sqrt(c);
  }
  else {
    double t=c-b*b/a;
    *pL11=sqrt(a);
    *pL21=b/ *pL11;
    
    if (t>0) {
      *pL22=sqrt(t);
    }
    else {
      assert(t > -0.000001); // numerical error
      *pL22=0;
    }
  }
}


void specsqrt(int n, double *pmx, double *psqrt) {
  gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(n);
  gsl_vector *eval = gsl_vector_alloc(n);  
  gsl_matrix *evec = gsl_matrix_alloc(n,n);
  double *ptmpmx=(double*) malloc(n*n*sizeof(double));
  gsl_matrix_view m = gsl_matrix_view_array(ptmpmx, n, n);
  int i,j;

// Otherwise MESSES UP pmx
  memcpy(ptmpmx,pmx,n*n*sizeof(double));

  gsl_eigen_symmv (&m.matrix, eval, evec, w);

  for (i=0;i<n;i++) {
    double *pthis= eval->data+i;
    if (*pthis<0) {
      //      printf("would have taken -ve sqrt of %f\n",*pthis);
      *pthis=0;
    }
    else {
      *pthis=sqrt(*pthis);
      //      printf("sqrt of eival is %f\n",*pthis);
    }
  }

  // GSL matrices are stored in row-major order (i.e. fill up 1st row first)
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      psqrt[i*n+j]=evec->data[i*n+j]*eval->data[j];
    }
  }

  free(ptmpmx);
  gsl_eigen_symmv_free(w);  
  gsl_vector_free(eval);
  gsl_matrix_free(evec);
}

double specdet(int n, double *pmx) {
  double *ptmpmx=(double*) malloc(n*n*sizeof(double));
  gsl_matrix_view m = gsl_matrix_view_array(ptmpmx, n, n);
  double eigenprod=1.0;
  int signum;
  gsl_permutation *p = gsl_permutation_alloc(n);

// Otherwise MESSES UP pmx
  memcpy(ptmpmx,pmx,n*n*sizeof(double));

  gsl_linalg_LU_decomp(&m.matrix,p,&signum);
  eigenprod = gsl_linalg_LU_det(&m.matrix, signum);

  free(ptmpmx);
  gsl_permutation_free(p);

  return(eigenprod);
}


void matrix_square_root(int n, double *pmx, double *psqrt) {
  if (n==1) {
    *psqrt=sqrt(*pmx);
  }
  else if (n==2) {
    chol22(pmx[0],pmx[2],pmx[3],&psqrt[0],&psqrt[2],&psqrt[3]);
    psqrt[1]=0;
  }
  else { // some eigenvalues could be zero so cannot do cholesky
    specsqrt(n, pmx, psqrt);
  }
}


double matrix_determinant(int n, double *pmx) {
  double det;
  if (n==1) {
    det=*pmx;
  }
  else if (n==2) {
    det=pmx[0] * pmx[3] - pmx[1]*pmx[2];
  }
  else {
    det=specdet(n,pmx);
  }
  if (det<0) {
    det=-det;
  }
  return(det);
}

int MaxOfIntVector(int n, int* pv) {
  int i,m=pv[0];

  for (i=1;i<n;i++) {
    if (pv[i]>m) {
      m=pv[i];
    }
  }
  return(m);
}

// Fills a matrix with the product of an n1xn2 and n3xn4, perhaps transposed
// C <- A x B

void mat_mult(double *pA, double *pB, double *pC, 
	      int n1, int n2, int n3,int n4, int T1, int T2) {
  int i,j,k;
  double sum;

  if (!T1 && !T2) {
    assert(n2==n3);
    for (i=0;i<n1;i++) {
      for (j=0;j<n4;j++) {
	sum=0.0;
	for (k=0;k<n2;k++) {
	  sum += pA[i*n2+k]*pB[k*n4+j];
	}
	pC[i*n4+j]=sum;
      }
    }
  }
  else if (!T1 && T2) {
    assert(n2==n4);
    for (i=0;i<n1;i++) {
      for (j=0;j<n3;j++) {
	sum=0.0;
	for (k=0;k<n2;k++) {
	  sum += pA[i*n2+k]* pB[j*n4+k];
	}
	pC[i*n3+j]=sum;
      }
    }
  }    
  else if (T1 && !T2) {
    assert(n1==n3);
    for (i=0;i<n2;i++) {
      for (j=0;j<n4;j++) {
	sum=0.0;
	for (k=0;k<n1;k++) {
	  sum += pA[k*n2+i]*pB[k*n4+j];
	}
	pC[i*n4+j]=sum;
      }
    }
  }
  else {
    assert(n1==n4);
    for (i=0;i<n2;i++) {
      for (j=0;j<n3;j++) {
	sum=0.0;
	for (k=0;k<n1;k++) {
	  sum += pA[k*n2+i]*pB[j*n4+k];
	}
	pC[i*n3+j]=sum;
      }
    }
  }
} 


void mat_mult_ld(double *pA, double *pB, double *pC, 
		 int n1, int n2, int n3,int n4, int T1, int T2) {
  int i,j,k;
  long double sum;

  if (!T1 && !T2) {
    assert(n2==n3);
    for (i=0;i<n1;i++) {
      for (j=0;j<n4;j++) {
	sum=0.0;
	for (k=0;k<n2;k++) {
	  sum += (long double)pA[i*n2+k]*(long double)pB[k*n4+j];
	}
	pC[i*n4+j]=(double)sum;
      }
    }
  }
  else if (!T1 && T2) {
    assert(n2==n4);
    for (i=0;i<n1;i++) {
      for (j=0;j<n3;j++) {
	sum=0.0;
	for (k=0;k<n2;k++) {
	  sum += (long double) pA[i*n2+k]*(long double) pB[j*n4+k];
	}
	pC[i*n3+j]=(double)sum;
      }
    }
  }    
  else if (T1 && !T2) {
    assert(n1==n3);
    for (i=0;i<n2;i++) {
      for (j=0;j<n4;j++) {
	sum=0.0;
	for (k=0;k<n1;k++) {
	  sum += (long double)pA[k*n2+i]*(long double)pB[k*n4+j];
	}
	pC[i*n4+j]=(double)sum;
      }
    }
  }
  else {
    assert(n1==n4);
    for (i=0;i<n2;i++) {
      for (j=0;j<n3;j++) {
	sum=0.0;
	for (k=0;k<n1;k++) {
	  sum += (long double)pA[k*n2+i]*(long double)pB[j*n4+k];
	}
	pC[i*n3+j]=(double)sum;
      }
    }
  }
} 


void vecvec_el_op(double *pA, double *pB, double *pC, int n, char op) {
  int i;
  switch(op) {
  case '+': {
    for (i=0;i<n;i++) {
      pC[i]=pA[i]+pB[i];
    }
    break;
  }
  case '-': {
    for (i=0;i<n;i++) {
      pC[i]=pA[i]-pB[i];
    }
    break;
  }
  case '*': {
    for (i=0;i<n;i++) {
      pC[i]=pA[i]*pB[i];
    }
    break;
  }
  case '/': {
    for (i=0;i<n;i++) {
      pC[i]=pA[i]/pB[i];
    }
    break;
  }
  default:
    printf("vecvec_el_op: unknown op %c\n",op);
    assert(7==8);
  }
}

double vecvec_dot_product(double *pA, double *pB, int n) {
  double *pC=(double*) malloc(n*sizeof(double));
  double res=0;
  int i;
  vecvec_el_op(pA, pB, pC, n, '*');
  for (i=0;i<n;i++) {
    res += pC[i];
  }
  free(pC);
  return(res);
}

void vec_el_op(double *pA, double *pC, int n, double val, char op) {
  int i;
  switch(op) {
  case '+': {
    for (i=0;i<n;i++) {
      pC[i]=pA[i]+val;
    }
    break;
  }
  case '-': {
    for (i=0;i<n;i++) {
      pC[i]=pA[i]-val;
    }
    break;
  }
  case '*': {
    for (i=0;i<n;i++) {
      pC[i]=pA[i]*val;
    }
    break;
  }
  case '/': {
    for (i=0;i<n;i++) {
      pC[i]=pA[i]/val;
    }
    break;
  }
  case 'e': {
    for (i=0;i<n;i++) {
      pC[i]=exp(pA[i]);
    }
    break;
  }
  case 'l': {
    for (i=0;i<n;i++) {
      pC[i]=log(pA[i]);
    }
    break;
  }
  default:
    printf("vec_el_op: unknown op %c\n",op);
    assert(7==8);
  }
}

void vec_set(double *pA, int n, double val) {
  int i;
  for (i=0;i<n;i++) {
    pA[i]=val;
  }
}


void mat_add_transpose(int n, double *pA, double *pC) {
  int i,j;
  double tmp;
  for (i=0;i<n;i++) {
    for (j=0;j<i;j++) {
      tmp=pA[n*i+j]+pA[n*j+i];
      pC[n*i+j]=tmp;
      pC[n*j+i]=tmp;
    }
    pC[n*i+i]=2*pA[n*i+i];
  }
}

void SymmetriseSquareMatrix(int n, double *px) {
  int i,j;
  double tmp;

  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      tmp=(px[i*n+j]+px[j*n+i])/2.0;
      px[i*n+j]=tmp;
      px[j*n+i]=tmp;
    }
  }

}


//int main(int argc, char *argv[]) {
int mainlinalg(void) {
  double pM1[9]={3,2,1,2,5,1,1,1,3};
  double pM2[25]={0.761300, 0.033724, 0.130690, 0.627355, -0.760266, 
-0.004834, 0.780219, -0.002107, -0.000347, 0.004908, 
0.126578, 0.025561, 4.668233, -2.012706, -0.126575, 
0.625926, 0.026753, -2.010806, 1.816321, -0.625917, 
-0.760270, -0.033651, -0.130690, -0.627349, 0.761285};
  //  double pM3[4]={1,2,3,4};
  double pM4[6]={5,6,7,8,9,10};
  double pM6[25];
  double pM7[2]={11,12};
  double pM8[3]={13,14,15};
  double pM9[20]={1,0,0,0,1, 0,1,0,0,0, 0,0,1,0,0, 0,0,0,1,0};
  double pI5[25]={1,0,0,0,0, 0,1,0,0,0, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,1};
  double pG[16]={0.205552, -0.019990, -0.098382, -0.205180, 
		 0.091748, 0.696232 ,-0.027709, -0.074057, 
		 -0.415112, 0.082863, 0.206637, 0.416911, 
		 -0.548763, 0.058057, 0.263802, 0.548110};
  double pPsi[16]={86266273327695488.000000, -1923872926233846.000000, -20666500820483464.000000, 96519517636481968.000000, 
		   -1923872926233846.000000, 42905377658321.796875, 460895317189423.000000, -2152536323534094.500000, 
		   -20666500820483464.000000, 460895317189423.000000, 4950999269152575.000000, -23122833681091224.000000, 
		   96519517636481968.000000, -2152536323534094.500000, -23122833681091224.000000, 107991419188806016.000000};
  double det;

  matrix_invert(3,pM1,pM6);
  printf("Inverted:\n");
  printdmat(3,3,pM1);
  printdmat(3,3,pM6);
  
  det=matrix_determinant(3,pM1);
  printf("Determinant %f\n",det);
  det=matrix_determinant(5,pM2);
  printf("Determinant %f\n",det);

  printf("multiplying 2x1 NN, 3x1 NT, 1x1TN, 1x2TT\n");
  printdmat(2,3,pM4);
  printdmat(1,2,pM7);
  printdmat(1,3,pM8);

  mat_mult(pM7,pM4,pM6,1,2,2,3,0,0);
  printdmat(1,3,pM6);  
  mat_mult(pM8,pM4,pM6,1,3,2,3,0,1);
  printdmat(1,2,pM6);
  mat_mult(pM4,pM4,pM6,2,3,2,3,1,0);
  printdmat(3,3,pM6);
  mat_mult(pM4,pM7,pM6,2,3,1,2,1,1);
  printdmat(3,1,pM6);

  printf("Pautoreg x identity\n");
  mat_mult(pM9,pI5,pM6,4,5,5,5,0,1);
  printdmat(4,5,pM6);

  printf("G Psi\n");
  mat_mult(pG,pPsi,pM6,4,4,4,4,0,0);
  printdmat(4,4,pG);
  printdmat(4,4,pPsi);
  printdmat(4,4,pM6);
  


  return(0);
}

int check_nonneg(double *pvec, int n) {
  int i=0, neg=0;

  while ((i<n) & !neg) {
    if (pvec[i]<0) {
      neg=1;
    }
    i++;
  }
  return(neg);
}
