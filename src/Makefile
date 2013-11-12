###########################################################
# Makefile
###########################################################
### Variables ###
#Created CSG
# Use GCC compiler
CC =  gfortran
CFLAGS =  -Wall -O3 -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF  
LDFLAGS =  -lgsl -lgslcblas -lm -lblas

base_dir=mcmc
sim_dir=$(base_dir)/simulators
#mcmc_dir=mcmc
generic_dir=$(base_dir)/generic

#LNA directories
codelib=$(sim_dir)/Generic_c
LinearAlgebra = $(codelib)/LinearAlgebra
Print = $(codelib)/Print
VGDIR=$(sim_dir)/VG_lsoda
mine=$(sim_dir)/LNA_generic_hybrid


vgdir_OBJECTS =	VG_bnorm.o VG_cfode.o VG_d1mach.o \
	VG_dgbfa.o VG_dgbsl.o VG_dgefa.o \
	VG_dgesl.o VG_ewset.o VG_fnorm.o \
	VG_hessian.o VG_intdy.o VG_integrate.o\
	VG_lsoda.o VG_prja.o VG_solsy.o\
	VG_stoda.o VG_vmnorm.o VG_xerrwx.o

mine_OBJECTS = struct_generic.o read.o system_spec.o \
	random.o hybrid.o genericODE.o 
linear_algebra_OBJECTS = cslinalg.o
print_OBJECTS = arrvec.o 
simulator_OBJECTS = hybrid_sde.o gillespie.o sde.o hybrid_lna.o
mcmc_OBJECTS = init.o likelihood.o Output.o update.o algorithm.o 
generic_OBJECTS=commandline.o timings.o paths.o

OBJECTS = $(foreach obj,$(mine_OBJECTS),$(mine)/$(obj))
OBJECTS += $(foreach obj,$(vgdir_OBJECTS),$(VGDIR)/$(obj))
OBJECTS += $(foreach obj,$(linear_algebra_OBJECTS),$(LinearAlgebra)/$(obj))
OBJECTS += $(foreach obj,$(print_OBJECTS),$(Print)/$(obj))
OBJECTS += $(foreach obj,$(simulator_OBJECTS),$(sim_dir)/$(obj)) 
OBJECTS += $(foreach obj,$(mcmc_OBJECTS),$(base_dir)/$(obj)) 
OBJECTS += $(foreach obj,$(generic_OBJECTS),$(generic_dir)/$(obj)) 

FILE=timing

#FILE=pmcmc
#pmcmc

DEBUG= #-pg 

$(FILE): $(FILE).c $(OBJECTS) 
	$(CC) $(CFLAGS) $(OBJECTS) -o $(FILE) $(FILE).c $(LDFLAGS) $(DEBUG)

%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $< $(DEBUG)

%.o : %.f
	$(CC) $(CFLAGS) -o $@ -c $< $(DEBUG)

clean:
	rm -fv  *~ $(FILE) $(OBJECTS) 
	rm -fv *.o

cleaner: clean
	rm -fv  output.csv parameters.csv residuals.csv prior.csv
