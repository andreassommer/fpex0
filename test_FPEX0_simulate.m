% simple tester for FPEX0_simulate

global FPEX0
FPEX0_setup()

% fraser suzuki parametrization
r  =   2;
h  =  40;
z  = 135.0;
wr =  15;
sr =   0.1;
p_IC = [r h z wr sr];

% FP drift and diffusion parametrization
p_FPdrift     = [ 0.00000  0.00000   ];
p_FPdiffusion = [ 0.00000  0.00000   ];

% assemble p
pvec = [ p_FPdrift  p_FPdiffusion  p_IC];

% simulate
sol = FPEX0_simulate(pvec);

% display solution
FPEX0_visualize(FPEX0, sol);