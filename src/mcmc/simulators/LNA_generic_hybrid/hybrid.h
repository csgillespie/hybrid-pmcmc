
void SetHybridInfo(HybridInfo *pHI, int *pA, double *px, double *ph, 
		   double *pc, double deltat);
void SetHybridRates(HybridInfo *pHI, double *ph);
void HybridModifyRates(double *ph, HybridInfo *pHI);
void HybridModifyDerivs(double *pdhdx, HybridInfo *pHI);

int *pISlow(HybridInfo *pHI);
double *pHybridSlowRates(HybridInfo *pHI);
int nSlow(HybridInfo *pHI);
double TotSlowRate(HybridInfo *pHI);

int IChooseSlowReaction(HybridInfo *pHI, double uniform);


