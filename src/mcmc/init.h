st_mcmc_settings * initMCMCSettings(int burn, int thin, long iters,
                                    char *dir, char *output_file);
st_model_at *initModelAttributes();
st_data * initData(st_model_at *model_at, char* dir);
st_part_at *initPart( st_parts_at *prior_parts, char *dir);
st_parts_at *initAcceptedParts(st_parts_at *prior_parts);
st_parts_at * initPrior(char *dir);

st_mcmc_update *initMCMCUpdate(char *dir);
