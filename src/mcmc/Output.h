void outputMCMCSettings(st_mcmc_settings *mcmc_settings);
void outputParameters(int iteration, st_mcmc_settings *mcmc_settings, st_part_at *accepted_part);
void outputData(st_data *data_at, int no_obs_sps);
void outputPartSettings(st_part_at *accepted_part);
void outputPartsSettings(st_parts_at *accepted_parts);
void outputParameterVector(st_part_at *part);
void outputTimings(double t_diff, st_mcmc_settings *mcmc_settings);
