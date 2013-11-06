void updateParameters(st_part_at *accepted_part, 
                      st_part_at *prop_parts, 
                      st_mcmc_update *mcmc_update);
void updateSpecies(int index, st_parts_at *prior_parts, st_part_at *prop_part);
void updateResiduals(int index, st_parts_at *prior_parts, st_part_at *prop_part);
void acceptParameters(st_model_at *model_at, st_part_at *accepted_part,
                      st_part_at *prop_part, 
                      int accept);
void switchParticles(st_parts_at *prior_parts, st_parts_at *accepted_parts);
void updateParticles(int iters, st_part_at *accepted_part, st_parts_at *accepted_parts);



