void get_A(int *pA, int S, int R);
void get_h(double *ph, double *px, double *pc, void *pparms);
void get_dhdx(double *pdhdx, double *px, double *pc, void *pparms);
void get_b(double *pb, double *px, double *pc);
void get_d2hdxmdxn(double *pd2hdxmdxn, double *px, double *pc, void *pparms,int hcpt);
void constrain(double *px);
void *get_fixed_parms();

void get_prior_for_X0(double *pmu, double *pSigma);

void initialise_rates(double *pc);
void initialise_obsvar(int Sobs, double *pvar);
int allowable_rates(double *pc);
int allowable_obsvar(int Sobs, double *pvar);
double log_prior_rates(double *pc);
double log_prior_obsvar(int Sobs, double *pobsvar);
void get_jumpvar(int Sobs, double *pjv);

void get_SR(int *pS, int *pR);
void get_nits_and_seed(int *pnits, int *pseed);
void get_nits_seed_thin(int *pnits, int *pseed, int *pthin);
char *pget_system();
char *pget_obsfname();
char *pget_outfname();
void get_hybrid_info(double *pN,double *phstar,double *pepsilon,
		     double *pdeltat, int *pdo_sk);

void free_fixed_parms(void *pparms);
