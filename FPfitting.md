# Approach

Reference:

\[1] Sommer, Hohenauer, Barz, 2022: Data-driven de-smearing of DSC data

## Optimization problem

See formula (4) in 

```
min_p0,pFP   sum_ij || y(t_j, x_i) -  fpeak(beta_j, T_i) ||^2

s.t.  d/dt  y(t,x) = -µ(t ; pFP) * d/dx y(t,x)  +  D(t ; pFP) d²/dx² y(t,x 

      y(0,x) = y0(x ; p0)      t \in [0, beta_max]   x \in IR

```

## Initial Distibution y0

Function of initial conditions should be freely selectable, but must be able to provide derivatives w.r.t. parameters p0.

So, we need:   `y0(x ; p0)`    and   `d/dp0 y0(x ; p0)`

Ideally, we have a certain collection of candidate functions, like scaled Gaussian, scaled Gompertz, scaled Weibull, Fraser-Suzuki, Haarhoff-Van-der-Linde, etc...

