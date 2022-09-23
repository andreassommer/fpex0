# FPEX0: Fokker-Planck-based extrapolation of DSC data to a zero heating rate


## Cloning the repository from github

This repository requires the following submodules (also available on github):
- mmtools
- optionlists

The easiest way to receive them is to use the `--recurse-submodules` option when cloning FPEX0:

```
git clone --recurse-submodules git@github.com:andreassommer/fpex0.git
```

Alternatively, you can download the required submodules manually from github into subfolders 
of the FPEX0 directory.




## About FPEX0

This repository gives a Matlab implementation of the FPEX0 method 
for data-driven de-smearing of DSC signals presented in the paper

Sommer, Andreas; Hohenauer, Wolfgang; Barz, Tilman:  
Data-driven de-smearing of DSC signals.  
J Therm Anal Calorim (2022).  
https://doi.org/10.1007/s10973-022-11258-y


## Running the test example.
0) Clone the repository, start Matlab and change to the folder where FPEX0 was downloaded to.
1)	Initialize necessary paths by invoking:   `FPEX0_initPaths();`
2)	If available, initialize a parallel cluster:   `pool=parpool(n);`  
    where `n` denotes the number of CPU cores available. 
    (Please see Matlab's documentation for details)
3)	Start the test problem by invoking: `FPEX0_exampleFit();`

The example uses Matlab's `lsqnonlin` for optimization; other optimizers are available (see `FPEX0_fit.m`). 

The example should converge within approx 20 iterations close to the following point:  
`x = [-0.9555 0.0328 0.2861 3.4172 2.5732 43.0439 131.8116 3.6591 0.1861]`

Depending on the problem and chosen accuracies, the optimization might stop due to different reasons.
A (close to) optimal solution is found, if the so called "first-order optimality" measure is small
(ideally zero, but when using FD approximations in FPEX0, values less than 10 are acceptable).



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
  
