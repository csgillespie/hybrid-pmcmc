/** Double **/

void printdvec(int n,double *pv);
void fprintdvec(int n,double *pv, FILE *pfp);
void fprintdvecnoret(int n,double *pv, FILE *pfp);
void printdmat(int nrow, int ncol,double *pv);
void fprintdmat(int nrow, int ncol, double *pv, FILE *pfp);
void printdarray(int nblock, int nrow, int ncol,double *pv);
void printdmat_precise(int nrow, int ncol,double *pv);

/** Integer **/

void printivec(int n,int *pv);
void fprintivec(int n, int *pv, FILE *pfp);
void printimat(int nrow, int ncol,int *pv);
void fprintimat(int nrow, int ncol, int *pv, FILE *pfp);
void printiarr(int nblock, int nrow, int ncol,int *pv);

void printstmat(int nrow, int ncol,size_t *pv);
