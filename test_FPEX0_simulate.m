% simple tester for FPEX0_simulate

FPEX0setup = FPEX0_exampleSetup();

% fraser suzuki parametrization
r  =     2;
h  =    40;
z  = 135.0;
wr =    15;
sr =   0.1;
p_IC = [r h z wr sr];

% FP drift and diffusion parametrization
p_FPdrift     = [ 0.00000  0.00000   ];
p_FPdiffusion = [ 0.00000  0.00000   ];

% assemble p
pvec = [ p_FPdrift  p_FPdiffusion  p_IC];

% NOTE: we could also use the initial parameter estimates from FPEX0setup
% pvec = FPEX0setup.p0_all;

% simulate
sol = FPEX0_simulate(FPEX0setup, pvec);

% display solution
FPEX0_visualize(FPEX0setup, sol);