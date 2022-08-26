# Approach

Reference:

\[1] Sommer, Hohenauer, Barz, 2022: Data-driven de-smearing of DSC data

## Optimization problem

See formula (4) in 

```
<<<<<<< HEAD
min_p0,pFP   sum_ij || y(t_j, x_i) -  f^peak(beta_j, T_i) ||^2

s.t.  d/dt  y(t,x) = -µ(t ; pFP) * d/dx y(t,x)  +  D(t ; pFP) d²/dx² y(t,x) 
=======
min_p0,pFP   sum_ij || y(t_j, x_i) -  fpeak(beta_j, T_i) ||^2

s.t.  d/dt  y(t,x) = -µ(t ; pFP) * d/dx y(t,x)  +  D(t ; pFP) d²/dx² y(t,x 
>>>>>>> main

      y(0,x) = y0(x ; p0)      t \in [0, beta_max]   x \in IR

```

<<<<<<< HEAD
where `y(t,x)` is the Fokker-Planck solution at time `t` and space `x`, and
`f^peak(beta_j, T_i)` denotes the measurements for heating rate `beta_j` and temperature `T_i`.

### Drift µ and diffusion D

Drift function `µ(t ; pFP)` and diffusion function `D(t ; pFP)` are linear functions.
With `pFP = (p^µ_1 , p^µ_2 , p^D_1 , p^D_2)`, we set:

```
µ(t ; pFP) = p^µ_1  +  p^µ_2 * t

```

To ensure a non-negative diffusion, a special formulation is used:

```
D(t ; pFP) = p^D_1  + t * (p^D_2 - p^D_1) / beta_max

```

By putting non-negativit constraints on `p^D_2` and `p^D_1`,
this ensures a non-negative diffusion in the interval `t \in [0, beta_max]`.

Note: A negative drift µ is valid and unproblematic!

### Initial Distibution y0
=======
## Initial Distibution y0
>>>>>>> main

Function of initial conditions should be freely selectable, but must be able to provide derivatives w.r.t. parameters p0.

So, we need:   `y0(x ; p0)`    and   `d/dp0 y0(x ; p0)`

Ideally, we have a certain collection of candidate functions, like scaled Gaussian, scaled Gompertz, scaled Weibull, Fraser-Suzuki, Haarhoff-Van-der-Linde, etc...

<<<<<<< HEAD
## Solving the PDE

The PDE is transformed into a (medium-scale) ODE by using the method of lines.
See FokkerPlanck1D.cpp in SolvIND for an example implementation.

We assume a Neumann boundary. That means that there might be some probability "loss"
at the borders, so the left and right boundaries must be chosen "far enough away" from
the probability support.

```
oneBYh = 1.0 / h;
oneBYhsquared = 1.0 / ( h*h );
dy[0] = oneBYhsquared * diffusion * (y[1] - y[0]);               // Neumann boundary (left)
for (j = 1; j < N-1; ++j) {
   dy[j] = - drift * oneBYh * 0.5 * (y[j+1] - y[j-1])                 // central difference
           + diffusion * oneBYhsquared * (y[j-1] - 2.0 * y[j] + y[j+1]);  // Laplace 3-star
   }
dy[N-1] = oneBYhsquared * diffusion * (y[N-2] - y[N-1]);        // Neumann boundary (right)

```

This results in an N-dimensional ODE, with initial value `y0 \in IR^N`

The ODE-solver should provide sensitivities w.r.t. the initial value, `dy/dy0`,
and w.r.t to the drift and diffusion parameters, `dy/dFP`. 

## Solving the optimization problem

For efficient solution, we need the derivatives w.r.t. the optimization variables
`pFP` and `p0`. We are only interested in the initial distribution `y0=y0(p0)`.
However, `pFP` must be adjusted for fitting.

So, what we actually need, are the following derivatives:
`dy/dpFP`  and  `dy/dp0 = dy/dy0  * dy0/dp0`.
`..(1)..`  and  `..(2)..= ..(3).. * ..(4)...`.

From the ODE solver, we get (1) and (3), i.e. sensitivities of the trajectory.

The sensitivity (4) `dy0/dp0` of the initial distribution is analytically computed
from the formula for the initial distriution.

###### Attention

Different initial distributions have a different number of parameters!
Also, drift and diffusion might be more complex parameterized.
So, the number of parameters must not be fixed to a specific value.
=======
>>>>>>> main
