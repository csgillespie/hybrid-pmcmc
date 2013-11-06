void part(st_model_at *model_at, 
          st_mcmc_settings *mcmc_settings, 
          void (*forwardSimulate)(st_part_at* prop_part, double tstep),
          char *dir);
