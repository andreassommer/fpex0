# FPEX0: Fokker-Planck-based extrapolation of DSC data to a zero heating rate

## Cloning the repository from github

This repository contains submodules.
The easiest way to receive them is to use the `--recurse-submodules` option:

```
git clone --recurse-submodules git@github.com:andreassommer/fpex0.git
```

## About FPEX0

This repository gives a Matlab implementation of the FPEX0 method 
for data-driven de-smearing of DSC signals presented in the paper

Sommer, Andreas; Hohenauer, Wolfgang; Barz, Tilman:  
Data-driven de-smearing of DSC signals.  
J Therm Anal Calorim (2022).  
https://doi.org/10.1007/s10973-022-11258-y


## Running the test example.
0)  Clone the repository, start Matlab and change to the folder where FPEX0 was downloaded to.
1)	Initialize necessary paths by invoking:   `FPEX0_initPaths();`
2)	If available, initialize a parallel cluster:   `pool=parpool(n);`  
    where `n` denotes the number of CPU cores available. 
    (Please see Matlab's documentation for details)
3)	Start the test problem by invoking: `FPEX0_exampleFit();`

The example uses Matlab's `lsqnonlin` for optimization; other optimizers are available (see `FPEX0_fit.m`). 

The example should converge to a solution within 29 iterations to the following solution:  
`x = [-0.9555 0.0328 0.2861 3.4172 2.6310 43.0438 131.8116 3.7368 0.1825]`



## Questions:
- Why is it slow?  
  Answer: Since the Matlab integrators are not capabale to deliver accurate and consistent sensitivity information,
  which is cruical for optimization, the required derivatives are approximated via external numerical
  differentiation. This requires on the one hand a lot of function evaluations, and on the other hand
  highly accurate integration tolerances, leading to lenghty computation times.
- Can it be made faster?  
  Answer: Yes. We have much faster external integrators at hand. Contact us for details.
- Further questions?  
  Contact us!  By email via andreas.sommer@iwr.uni-heidelberg.de
  
