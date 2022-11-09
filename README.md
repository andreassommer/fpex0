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
`p = [-0.9553 0.0329 0.2693 3.4700 2.5736 42.0826 131.7993 3.7684 0.1909]`

Depending on the problem and chosen accuracies, the optimization might stop due to different reasons.
A (close to) optimal solution is found, if the so called "first-order optimality" measure is small
(ideally zero, but when using FD approximations in FPEX0, values less than 10 are acceptable. Default is 1.0).



## Questions:
- Why is it slow?  
  Answer: Since the Matlab integrators are not capabale to deliver accurate and consistent sensitivity information,
  which is cruical for optimization, the required derivatives are approximated via (1) external numerical differentiation,
  or by using the (2) variational differential equations.  
  (1) requires lot of function evaluations, and highly accurate integration tolerances, leading to lenghty computation times. 
  (2) is much faster, and in principle more accurate, but does not deliver consistent sensitivities.
  With state-of-the-art solvers like SolvIND from our group, we can deliver consistent and highly accurate derivatives.
- Can it be made faster?  
  Answer: Yes. We have much faster external integrators at hand. Contact us for details.
- Can I do anything to make it faster?  
  There are several possibilities:  
  - Reduce the integration tolerances, e.g. set reltol = 1e-6 and abstol = 1-e10. 
    When using the VDE sensitivities (default), the tolerances for the nominal might even be lower.
  - Change finite difference approximation in the optimizer settings to "forward"  
  - Use a coarser temperature grid - or a finer! The coarser, the more stiff the system becomes at initial times.  
  Note that the aforementioned tricks might negatively impact the optimizer's optimality check, 
  making it hard to detect a "mathematically clean" optimum by testing the first-order optimality,
  but should still converge to a "good" fit.  
  This fit can then be re-evaluated with finer grids, tighter tolerances, etc., to show optimality (if required).
- Further questions?  
  Contact us!  By email via andreas.sommer@iwr.uni-heidelberg.de or code@andreas-sommer.eu.


