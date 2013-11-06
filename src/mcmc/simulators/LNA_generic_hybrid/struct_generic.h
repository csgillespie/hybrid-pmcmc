typedef struct _reaction_info {
  int S; // number of species
  int R; // number of reactions
  int maxnreact; //  max number of reactions a species can be involved in
  int maxnspec; // max number of species affected net by a reaction
  double *pc; // rate constants
  int *pA; // Net effects matrix
  int *pNS; // number of species AFFECTED (net) by a given reaction
  int *pS;  // species affected (net) by a given reaction
  int *pNR; // # reactions that change the net amount, for each species
  int *pR; // reactions which change the net amount, for each species
  int *pNRR; // # reactions which change the net amount of species i and j
  int *pRR; // reactions which change the net amount of species i and j
  void *pparms;
} ReactionInfo;

typedef struct _hybrid_info {
  int R; // number of reactions
  int S; // number of species
  int *pifast; // int indicator 1=fast
  int *pislow; // int indicator 1=slow
  int nslow;
  int nfast;
  int hybridrun;
  double *pslow_rates; // rate for slow but zero for fast
  double *pslow_c; // c for slow but zero for fast
  double tot_slow;
  double max_lamslow;
  double *pbstar;
  double *pbstarmax; // running maximum of bstars
  double N; // S&K parameter - need x_{ij} > N |a_{ij}| for fast
  double hstar; // S&K parameter - need h_j > hstar Deltat for fast
  double epsilon; // us parameter
  double deltat; // for hybrid decisions
  int do_sk;
} HybridInfo;


typedef struct _sde_info {
  int S;
  int R;
  double *ph;
  double *pmu;
  double *pF;
  double *pS2;
  double *pS;
  double *pHaV;
  double *pHS2V;
  double *pscratch;
  int keepmaxs;
  ReactionInfo *pRI;
  HybridInfo *pHI;
} SDEInfo;


typedef struct _lna_realise {
  int S;
  double t;
  double *py; // ODE
  double *pM; // discrepancy from deterministic value 
  double *px; // realisation
} LNARealisation;

typedef struct _obs {
  int S;
  int Sobs;
  double t;
  double *pobs; 
  double *pP; // converts: pobs = P X
} Observation;

typedef struct _meanvar {
  int S;
  double t; 
  double *pmu;
  double *pSigma;
} MeanVar; // mean and variance for M_t | obs up to t-1

// ** Print **

void PrintReactionInfo(ReactionInfo *pI);
void PrintSDEInfo(SDEInfo *pSI);
void PrintLNARealisation(LNARealisation *pLR);
void PrintObservation(Observation *pOI);
void PrintMeanVar(MeanVar *pMV);

void PrintLNARealisationArray(int n, LNARealisation **ppRI);
void PrintObservationArray(int n, Observation **ppOI);
void PrintMeanVarArray(int n, MeanVar **ppMV);

void FPrintRealisation(LNARealisation *pLR, FILE *pfp);
void FPrintAllRealisations(int n, LNARealisation **ppRI, FILE *pfp);

void PrintHybridInfo(HybridInfo *pHI);



// ********** Allocate and Free ************

void AllocSDEInfo(ReactionInfo *pI, SDEInfo *pSI);
void DeAllocSDEInfo(SDEInfo *pSI);
void AllocLNARealisation(ReactionInfo *pRI, LNARealisation *pLR);
void DeAllocLNARealisation(LNARealisation *pLR);
void AllocObservation(int S, int Sobs, Observation *pOI, double *pP);
void DeAllocObservation(Observation *pOI);
void AllocMeanVar(ReactionInfo *pRI, MeanVar *pMV);
void DeAllocMeanVar(MeanVar *pMV);

void AllocLNARealisationArray(int n_times, ReactionInfo *pRI, 
			      LNARealisation **ppLR);
void AllocObservationArray(int n_times, int S, int Sobs, Observation **ppOI,
			   double *pP);
void AllocMeanVarArray(int n_times, ReactionInfo *pRI, MeanVar **ppMV);
void DeAllocLNARealisationArray(int n_times, LNARealisation **ppLR);
void DeAllocObservationArray(int n_times, Observation **ppOI);
void DeAllocMeanVarArray(int n_times, MeanVar **ppMV);

void AllocHybridInfo(ReactionInfo *pI, HybridInfo *pHI);
void DeAllocHybridInfo(HybridInfo *pHI);


