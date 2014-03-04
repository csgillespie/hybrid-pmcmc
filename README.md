Bayesian Inference for Hybrid Discrete-Continuous Systems Biology Models
------------------------------------------------------------------------
  
This repository contains the source code for the paper

[Sherlock, C.](http://www.maths.lancs.ac.uk/~sherlocc/), [Golightly, A.](http://www.mas.ncl.ac.uk/~nag48/), [Gillespie, C.S.](http://www.mas.ncl.ac.uk/~ncsg3/) (2013) *Bayesian Inference for Hybrid Discrete-Continuous Systems Biology Models* 
([arXiv](http://arxiv.org/abs/1402.6602))


To download the code, either clone this repository or download the [zip](https://github.com/csgillespie/hybrid-pmcmc/archive/master.zip) file and unpack.

1. Change to the `src` directory
```
cd src
```

1. Type 
```
make
```
the code should build with no errors and create an execuable called `pmcmc`. Note the code the `gsl` math library is needed.

1. The `pmcmc` has a number of options:
 * `-s` the simulator used - either `gillespie`, `hybridLNA`, `hybridSDE`
 * `-n` the number of iterations
 * `-d` the input directory. This directory should contain a sub-directory with the name of the simulator used.
 * `-b` the burnin - default 0
 * `-t` thin - default 1

  For example
```
./pmcmc -s gillespie -n 1000 -d ../input/inference/run1/ -t 2
```

1. Running the code will create two files:

 * `name_of_simulator.csv` - mcmc output;
 * `timing.csv` - run time (in seconds)


