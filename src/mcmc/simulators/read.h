void ReadEffectsMatrices(char *pfname, int **ppL, int **ppR, int *pS, int *pR);
void ReadNObsAndSobs(FILE *pfp, int *pnobs, int *pSobs);
void ReadAllObs(FILE *pfp, int maxlen, int nobs, Observation **ppOI);
void ReadP(FILE *pfp, int maxlen, int S, int Sobs, double *pP);
