void matrix_invert(int n, double *pmx, double *pinv); // square
void matrix_square_root(int n, double *pmx, double *psqrt); // square
double matrix_determinant(int n, double *pmx); //square
int MaxOfIntVector(int n, int* pv);
void mat_mult(double *pA, double *pB, double *pC, int n1, int n2, int n3,int n4,
	      int T1, int T2);
void vecvec_el_op(double *pA, double *pB, double *pC, int n, char op);
double vecvec_dot_product(double *pA, double *pB, int n);
void vec_el_op(double *pA, double *pC, int n, double val, char op);
void vec_set(double *pA, int n, double val);
void mat_add_transpose(int n, double *pA, double *pC); // square
void SymmetriseSquareMatrix(int n, double *px);
