% simple tester for FPEX0_simulate

FPEX0setup = FPEX0_exampleSetup();

% integrator options
FPEX0setup.Integration.options.Stats = 'on';

% fraser suzuki parametrization
r  =     2;
h  =    40;
z  = 135.0;
wr =    15;
sr =   0.1;
p_IC = [r h z wr sr];

% FP drift and diffusion parametrization
p_FPdrift     = [  1.0      -0.0001   ];
p_FPdiffusion = [  1.001   0.001   ];

% assemble p
pvec = [ p_FPdrift  p_FPdiffusion  p_IC];

% simulate
tic
sol = FPEX0_simulate(FPEX0setup, pvec);
toc

% display solution
FPEX0_visualize(FPEX0setup, sol);