source("stoc_simulation.R")
load("choleskis.RData")
##Globals
estimate_particles = TRUE
RUN = 1
simulator = "gillespie"


seed = 3
sc = 10^(0:4)[RUN]
pars =  c(2,sc,1/50,1,1/(50*sc))

##Estimate equilibrium values
f = function(x) {
  x1 = x[1];x2=x[2]
  ((2-x1/50 - x1*x2 /(50*sc))^2 + (sc-x2+ 20*x1*x2/(50*sc))^2)
}


ic = floor(nlm(f, c(10, 10))$estimate/2)
ic = c(0, 0)
parts_matrix = matrix(c(300, 300, 1750,
                        800, 800, 1500, 
                        65, 65, 125, 
                        65, 65, 85,
                        30, 30, 40), ncol=3, byrow=T)
colnames(parts_matrix) = c("gillespie", "hybridLNA", "hybridSDE")
no_parts = parts_matrix[RUN,simulator]

#dir =  paste0("../input/run", RUN, "/")
if(estimate_particles) {
  input_dir =  paste0("../input/estimate_particles/run", RUN, "/", simulator, "/")
}else{
  input_dir =  paste0("../input/inference/run", RUN, "/", simulator, "/")
}

## R code for constructing files for the SIR model
set.seed(seed)

#Species Prior
df1 = data.frame(rep(ic[1], no_parts), ic[2])
write.table(df1, paste0(input_dir, "species_prior.csv"), sep=",", row.names=F, col.names=F)


#Parameter Prior
p1 = pars[1]#exp(runif(1, -10, 3))
p2 = pars[2]#exp(runif(1, -10, 3))
p3 = pars[3]#exp(runif(no_parts, -5, 3))
p4 = pars[4]#exp(runif(1, -10, 3))
p5 = pars[5]#exp(runif(1, -10, 3))
df2 = round(data.frame(p1, p2, p3, p4, p5), 6)
write.table(df2, paste0(input_dir, "parameters_prior.csv"), sep=",", row.names=F, col.names=F)

#Residual Prior
p1 = log(runif(no_parts))
p2 = log(runif(no_parts))
p3 = log(runif(no_parts))
p4 = log(runif(no_parts))
p5 = log(runif(no_parts))
df3 = data.frame(p1, p2, p3, p4, p5)
write.table(df3, paste0(input_dir, "residuals.csv"), sep=",", row.names=F, col.names=F)


write.table(choleskis[[RUN]], 
            file=paste0(input_dir, "tuning.csv"), sep=",", 
            col.names=F, row.names=F)

## Fix parameters
fix = matrix(0, ncol=5); fix[3] = 1;
if(estimate_particles) fix = fix + 1
#fix[5] = 1

write.table(fix, file=paste0(input_dir, "fixed.csv"), sep=",", col.names=F, row.names=F)


l = generate_simulation(ic, pars, 1, 50, seed=seed, filename=paste0(input_dir, "full.csv"))


#################################################
#Save settings
#################################################
settings = paste0(input_dir, "settings.R")
system(paste("cp input.R", settings))

